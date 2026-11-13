"""The window."""

from __future__ import annotations

import os
import queue
import re
import subprocess
import sys
import threading
import time
import tkinter as tk
import webbrowser
import zipfile
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from varan_gui import __version__, jobs
from varan_gui.engine import DOWNLOAD_PAGE, Engine
from varan_gui.settings import Settings, app_dir, cancer_types, resource, vep_cache_versions

ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
TSV = [("Tables", "*.tsv *.txt *.csv"), ("All files", "*.*")]

# Log lines that start a new step of the run (varan.py, walk.py, the
# update/extract/remove scripts) and what the status bar says for them.
STAGES = [
    ("Starting preparation study folder", "Preparing the study folder"),
    ("Starting vcf2maf conversion of", "Annotating variants with VEP"),
    ("Starting MAF filtering", "Filtering variants"),
    ("Concatenating mutation file", "Merging the mutations"),
    ("It's time to create tables!", "Writing the cBioPortal tables"),
    ("Starting validation", "Validating the study"),
    ("Starting Update study", "Merging the studies"),
    ("Starting Extract sample", "Extracting the samples"),
    ("Starting Delete sample", "Removing the samples"),
]


def open_path(path: str | Path) -> None:
    if sys.platform == "win32":
        os.startfile(str(path))  # noqa: S606
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


class PathField:
    """Label, entry and "Browse..." button on one grid row."""

    def __init__(self, parent, row: int, label: str, var: tk.StringVar, kind: str = "file",
                 filetypes=None, app: App | None = None, optional: bool = False) -> None:
        text = f"{label} (optional)" if optional else label
        ttk.Label(parent, text=text).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        self.entry = ttk.Entry(parent, textvariable=var, width=60)
        self.entry.grid(row=row, column=1, sticky="ew", pady=3)
        self.button = ttk.Button(parent, text="Browse...", command=self.browse)
        self.button.grid(row=row, column=2, padx=(6, 0), pady=3)
        self.var, self.kind, self.filetypes, self.app = var, kind, filetypes or [("All files", "*.*")], app

    def browse(self) -> None:
        start = self.var.get().strip() or (self.app.settings.last_dir if self.app else "")
        start_dir = start if Path(start).is_dir() else str(Path(start).parent) if start else ""
        if self.kind == "dir":
            value = filedialog.askdirectory(initialdir=start_dir or None, mustexist=True)
        else:
            value = filedialog.askopenfilename(initialdir=start_dir or None, filetypes=self.filetypes)
        if value:
            self.var.set(str(Path(value)))
            if self.app:
                self.app.settings.last_dir = str(Path(value) if self.kind == "dir" else Path(value).parent)

    def enable(self, on: bool) -> None:
        state = "normal" if on else "disabled"
        self.entry.configure(state=state)
        self.button.configure(state=state)


class Scrollable(ttk.Frame):
    """A frame whose content scrolls vertically."""

    def __init__(self, parent) -> None:
        super().__init__(parent)
        canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0)
        bar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.body = ttk.Frame(canvas, padding=12)
        self.body.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        window = canvas.create_window((0, 0), window=self.body, anchor="nw")
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        canvas.configure(yscrollcommand=bar.set)
        canvas.pack(side="left", fill="both", expand=True)
        bar.pack(side="right", fill="y")
        self.canvas = canvas
        self.body.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", self._wheel))
        self.body.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))

    def _wheel(self, event) -> None:
        step = -1 if event.delta > 0 else 1
        self.canvas.yview_scroll(step * (1 if sys.platform == "darwin" else 3), "units")


class Intro(tk.Toplevel):
    """The intro video (assets/intro_frames.zip, see make_intro.py) in a
    borderless window; a click or Esc skips it."""

    def __init__(self, master, frames: Path, on_done) -> None:
        super().__init__(master)
        self.archive = zipfile.ZipFile(frames)
        self.names = sorted(n for n in self.archive.namelist() if n.endswith(".gif"))
        try:
            self.fps = float(self.archive.read("fps.txt"))
        except (KeyError, ValueError):
            self.fps = 12.0
        self.on_done = on_done
        self.overrideredirect(True)
        self.configure(background="black")
        self.attributes("-topmost", True)
        self.screen = tk.Label(self, borderwidth=0, background="black", cursor="hand2")
        self.screen.pack()
        # tk images are in screen pixels: on a 150% / 200% display the frames
        # are enlarged (zoom, then subsample) to look as big as on a 100% one.
        scale = master.winfo_fpixels("1i") / 96
        self.zoom = (1, 1) if scale < 1.25 else (3, 2) if scale < 1.75 else (2, 1)
        self.frame = self._load(self.names[0])
        self.screen.configure(image=self.frame)
        w, h = self.frame.width(), self.frame.height()
        self.geometry(f"{w}x{h}+{(self.winfo_screenwidth() - w) // 2}+{(self.winfo_screenheight() - h) // 2}")
        for widget in (self, self.screen):
            widget.bind("<Button-1>", lambda e: self.finish())
        self.bind("<Escape>", lambda e: self.finish())
        self.focus_force()
        self.started = time.monotonic()
        self.shown = 0
        self.after(int(1000 / self.fps), self._next)

    def _next(self) -> None:
        if not self.winfo_exists():
            return
        # Paced on the clock: a slow computer drops frames instead of
        # stretching the video.
        index = int((time.monotonic() - self.started) * self.fps)
        if index >= len(self.names):
            self.finish()
            return
        if index != self.shown:
            self.frame = self._load(self.names[index])
            self.screen.configure(image=self.frame)
            self.shown = index
        self.after(max(1, int(1000 / self.fps / 2)), self._next)

    def _load(self, name: str) -> tk.PhotoImage:
        image = tk.PhotoImage(data=self.archive.read(name))
        zoom, sub = self.zoom
        if zoom > 1:
            image = image.zoom(zoom)
        return image.subsample(sub) if sub > 1 else image

    def finish(self) -> None:
        if self.on_done is None:
            return
        done, self.on_done = self.on_done, None
        self.archive.close()
        self.destroy()
        done()


def section(parent, title: str, row: int) -> ttk.LabelFrame:
    frame = ttk.LabelFrame(parent, text=title, padding=10)
    frame.grid(row=row, column=0, sticky="ew", pady=(0, 10))
    frame.columnconfigure(1, weight=1)
    return frame


class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.settings = Settings.load()
        self.events: queue.Queue = queue.Queue()
        self.engine = Engine(self.settings.engine)
        self.engine_ready = False
        self.busy = False
        self.proc: subprocess.Popen | None = None
        self.run: jobs.Run | None = None
        self.result: Path | None = None
        self.stopping = False

        self.title("Varan")
        # Sizes in screen pixels, scaled like the text on 150% / 200% displays.
        scale = max(1.0, self.winfo_fpixels("1i") / 96)
        width = min(int(1040 * scale), self.winfo_screenwidth() - 40)
        height = min(int(860 * scale), self.winfo_screenheight() - int(90 * scale))
        self.geometry(f"{width}x{height}+{(self.winfo_screenwidth() - width) // 2}+{int(20 * scale)}")
        self.minsize(int(820 * scale), int(600 * scale))
        try:
            self.icon = tk.PhotoImage(file=str(resource("assets/varan.png")))
            self.iconphoto(True, self.icon)
        except tk.TclError:
            self.icon = None

        # The intro plays while the window is built and Docker is checked.
        frames = resource("assets/intro_frames.zip")
        if self.settings.intro == "yes" and frames.is_file():
            self.withdraw()
            try:
                Intro(self, frames, on_done=self._show).update()
            except (OSError, zipfile.BadZipFile, tk.TclError, IndexError):
                self.deiconify()

        self._build()
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.after(100, self._poll)
        self.check_engine(start=True)

    def _show(self) -> None:
        self.deiconify()
        self.lift()
        self.focus_force()

    # ------------------------------------------------------------- layout

    def _build(self) -> None:
        style = ttk.Style(self)
        if "vista" not in style.theme_names() and "aqua" not in style.theme_names():
            style.theme_use("clam")
        style.configure("Title.TLabel", font=("TkDefaultFont", 16, "bold"))
        style.configure("Run.TButton", padding=(14, 6))

        head = ttk.Frame(self, padding=(12, 10, 12, 4))
        head.pack(fill="x")
        try:
            self.logo = tk.PhotoImage(file=str(resource("assets/varan_header.png")))
            scale = self.winfo_fpixels("1i") / 96
            if scale >= 1.25:
                self.logo = self.logo.zoom(3).subsample(2) if scale < 1.75 else self.logo.zoom(2)
            ttk.Label(head, image=self.logo).pack(side="left", padx=(0, 10))
        except tk.TclError:
            self.logo = None
        titles = ttk.Frame(head)
        titles.pack(side="left")
        ttk.Label(titles, text="Varan", style="Title.TLabel").pack(anchor="w")
        ttk.Label(titles, text="cBioPortal studies from sequencing results").pack(anchor="w")
        engine_box = ttk.Frame(head)
        engine_box.pack(side="right")
        self.engine_label = ttk.Label(engine_box, text="Checking Docker...")
        self.engine_label.pack(side="left", padx=(0, 8))
        self.engine_button = ttk.Button(engine_box, text="Retry", command=lambda: self.check_engine(start=True))
        self.engine_button.pack(side="left")
        self.engine_button.pack_forget()

        panes = ttk.PanedWindow(self, orient="vertical")
        panes.pack(fill="both", expand=True, padx=12, pady=(4, 0))
        self.tabs = ttk.Notebook(panes)
        panes.add(self.tabs, weight=3)
        self._build_create()
        self._build_update()
        self._build_samples("extract", "Extract samples")
        self._build_samples("remove", "Remove samples")
        self._build_settings()
        panes.add(self._build_log(panes), weight=2)
        self.after(300, lambda: panes.sashpos(0, int(panes.winfo_height() * 0.68)))

    def _buttons(self, parent, row: int, run, check=None) -> None:
        bar = ttk.Frame(parent)
        bar.grid(row=row, column=0, sticky="e", pady=(4, 0))
        if check:
            ttk.Button(bar, text="Check inputs", command=check).pack(side="left", padx=(0, 8))
        button = ttk.Button(bar, text="Run Varan", style="Run.TButton", command=run)
        button.pack(side="left")
        self.run_buttons = [*getattr(self, "run_buttons", []), button]

    def _output_fields(self, frame, row: int, prefix: str) -> None:
        out_dir = tk.StringVar(value=self.settings.last_dir)
        name = tk.StringVar()
        setattr(self, f"{prefix}_out_dir", out_dir)
        setattr(self, f"{prefix}_name", name)
        PathField(frame, row, "Output folder", out_dir, kind="dir", app=self)
        ttk.Label(frame, text="Study folder name").grid(row=row + 1, column=0, sticky="w", padx=(0, 12), pady=3)
        ttk.Entry(frame, textvariable=name, width=30).grid(row=row + 1, column=1, sticky="w", pady=3)
        ttk.Label(frame, text='Varan adds the version: "name_v1", "name_v2"...',
                  foreground="gray40").grid(row=row + 2, column=1, sticky="w")

    def _build_create(self) -> None:
        page = Scrollable(self.tabs)
        self.tabs.add(page, text="New study")
        body = page.body
        body.columnconfigure(0, weight=1)
        v = self.create_vars = {
            "input_mode": tk.StringVar(value="table"),
            "sample_tsv": tk.StringVar(), "patient_tsv": tk.StringVar(), "fusion_tsv": tk.StringVar(),
            "input_folder": tk.StringVar(), "pipeline": tk.StringVar(value=jobs.PIPELINES[""]),
            "sample_type": tk.StringVar(value="Solid"), "cancer": tk.StringVar(),
            "project_id": tk.StringVar(), "project_name": tk.StringVar(), "description": tk.StringVar(),
            "existing": tk.StringVar(value="new"), "oncokb": tk.BooleanVar(),
            "multiple": tk.BooleanVar(), "sigma": tk.BooleanVar(),
        }
        self.filter_vars = {k: tk.BooleanVar() for k in jobs.FILTERS}

        box = section(body, "Input", 0)
        modes = ttk.Frame(box)
        modes.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        ttk.Radiobutton(modes, text="Sample table (sample.tsv)", value="table",
                        variable=v["input_mode"], command=self._input_mode).pack(side="left", padx=(0, 16))
        ttk.Radiobutton(modes, text="Input folder", value="folder",
                        variable=v["input_mode"], command=self._input_mode).pack(side="left")
        self.table_fields = [
            PathField(box, 1, "Sample table", v["sample_tsv"], filetypes=TSV, app=self),
            PathField(box, 2, "Patient table", v["patient_tsv"], filetypes=TSV, app=self, optional=True),
            PathField(box, 3, "Fusion table", v["fusion_tsv"], filetypes=TSV, app=self, optional=True),
        ]
        self.folder_field = PathField(box, 4, "Input folder", v["input_folder"], kind="dir", app=self)
        templates = ttk.Label(box, text="Templates of the tables", foreground="#1a6fb3", cursor="hand2")
        templates.grid(row=5, column=1, sticky="w")
        templates.bind("<Button-1>", lambda e: webbrowser.open(
            "https://github.com/bioinformatics-policlinicogemelli/Varan/tree/main/Templates"))

        box = section(body, "Assay", 1)
        ttk.Label(box, text="Pipeline").grid(row=0, column=0, sticky="w", padx=(0, 12), pady=3)
        pipe = ttk.Combobox(box, textvariable=v["pipeline"], state="readonly", width=40,
                            values=list(jobs.PIPELINES.values()))
        pipe.grid(row=0, column=1, sticky="w", pady=3)
        pipe.bind("<<ComboboxSelected>>", lambda e: self._pipeline())
        ttk.Label(box, text="Sample type").grid(row=1, column=0, sticky="w", padx=(0, 12), pady=3)
        self.type_box = ttk.Combobox(box, textvariable=v["sample_type"], state="readonly", width=12,
                                     values=["Solid", "Liquid"])
        self.type_box.grid(row=1, column=1, sticky="w", pady=3)

        box = section(body, "Study", 2)
        ttk.Label(box, text="Cancer type").grid(row=0, column=0, sticky="w", padx=(0, 12), pady=3)
        self.cancers = [f"{cid} - {name}" for cid, name in cancer_types()]
        cancer = ttk.Combobox(box, textvariable=v["cancer"], values=self.cancers, width=45)
        cancer.grid(row=0, column=1, sticky="w", pady=3)
        cancer.bind("<KeyRelease>", lambda e: self._filter_cancers(cancer))
        ttk.Label(box, text='Type to search, e.g. "lung" - or "mixed" for several tumour types',
                  foreground="gray40").grid(row=1, column=1, sticky="w")
        for row, (key, label) in enumerate((("project_id", "Study id"), ("project_name", "Study name"),
                                            ("description", "Description")), start=2):
            ttk.Label(box, text=f"{label} (optional)").grid(row=row, column=0, sticky="w", padx=(0, 12), pady=3)
            ttk.Entry(box, textvariable=v[key], width=60).grid(row=row, column=1, sticky="ew", pady=3)

        box = section(body, "Output", 3)
        self._output_fields(box, 0, "create")
        ttk.Label(box, text="If the study exists").grid(row=3, column=0, sticky="nw", pady=(6, 3))
        modes = ttk.Frame(box)
        modes.grid(row=3, column=1, sticky="w", pady=(6, 3))
        for value, label in (("new", "make a new version"), ("overwrite", "overwrite the last version"),
                             ("resume", "resume the last version (after a stop or an error)")):
            ttk.Radiobutton(modes, text=label, value=value, variable=v["existing"]).pack(anchor="w")

        box = section(body, "Analysis", 4)
        ttk.Checkbutton(box, text="OncoKB annotation (needs a token in Settings)",
                        variable=v["oncokb"]).grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Checkbutton(box, text="Multi-sample VCF files", variable=v["multiple"]).grid(
            row=1, column=0, columnspan=3, sticky="w")
        ttk.Checkbutton(box, text="SigMA mutational signatures (HRD, Signature 3)",
                        variable=v["sigma"]).grid(row=2, column=0, columnspan=3, sticky="w")
        ttk.Label(box, text="Variant filters").grid(row=3, column=0, sticky="nw", pady=(8, 0))
        grid = ttk.Frame(box)
        grid.grid(row=3, column=1, columnspan=2, sticky="w", pady=(8, 0))
        for n, (letter, label) in enumerate(jobs.FILTERS.items()):
            ttk.Checkbutton(grid, text=label, variable=self.filter_vars[letter]).grid(
                row=n // 2, column=n % 2, sticky="w", padx=(0, 18))
        ttk.Label(box, text="Thresholds and lists come from conf.ini (Settings).",
                  foreground="gray40").grid(row=4, column=1, sticky="w", pady=(4, 0))

        self._buttons(body, 5, run=lambda: self.start(self._create_form()),
                      check=lambda: self.start(self._create_form(), check=True))
        self._input_mode()

    def _build_update(self) -> None:
        page = Scrollable(self.tabs)
        self.tabs.add(page, text="Merge studies")
        body = page.body
        body.columnconfigure(0, weight=1)
        self.update_study, self.update_new = tk.StringVar(), tk.StringVar()
        self.update_overwrite = tk.BooleanVar()
        box = section(body, "Studies", 0)
        PathField(box, 0, "Study", self.update_study, kind="dir", app=self)
        PathField(box, 1, "Study to add", self.update_new, kind="dir", app=self)
        ttk.Label(box, text="Samples of the second study are added to the first; a sample in both "
                  "is taken from the second.", foreground="gray40").grid(row=2, column=1, sticky="w")
        box = section(body, "Output", 1)
        self._output_fields(box, 0, "update")
        ttk.Checkbutton(box, text="Overwrite the last version", variable=self.update_overwrite).grid(
            row=3, column=1, sticky="w", pady=(6, 0))
        self._buttons(body, 2, run=lambda: self.start(self._study_form("update")))

    def _build_samples(self, operation: str, title: str) -> None:
        page = Scrollable(self.tabs)
        self.tabs.add(page, text=title)
        body = page.body
        body.columnconfigure(0, weight=1)
        study, overwrite = tk.StringVar(), tk.BooleanVar()
        box = section(body, "Study", 0)
        PathField(box, 0, "Study", study, kind="dir", app=self)
        box = section(body, "Samples to " + operation, 1)
        listbox = tk.Listbox(box, selectmode="extended", height=10, exportselection=False)
        listbox.grid(row=0, column=0, columnspan=2, sticky="ew")
        bar = ttk.Scrollbar(box, orient="vertical", command=listbox.yview)
        bar.grid(row=0, column=2, sticky="ns")
        listbox.configure(yscrollcommand=bar.set)
        box.columnconfigure(0, weight=1)
        tools = ttk.Frame(box)
        tools.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 0))
        count = ttk.Label(tools, text="Choose a study to list its samples.", foreground="gray40")
        ttk.Button(tools, text="Select all", command=lambda: (listbox.selection_set(0, "end"), refresh())).pack(
            side="left", padx=(0, 6))
        ttk.Button(tools, text="Clear", command=lambda: (listbox.selection_clear(0, "end"), refresh())).pack(
            side="left", padx=(0, 6))
        ttk.Button(tools, text="Select from a list file...", command=lambda: from_file()).pack(
            side="left", padx=(0, 12))
        count.pack(side="left")

        def refresh(*_) -> None:
            n = len(listbox.curselection())
            total = listbox.size()
            count.configure(text=f"{n} of {total} selected" if total else "No samples found in this folder.")

        def load(*_) -> None:
            listbox.delete(0, "end")
            folder = Path(study.get().strip()) if study.get().strip() else None
            if folder is None or not folder.is_dir():
                count.configure(text="Choose a study to list its samples.")
                return
            for sample in jobs.study_samples(folder):
                listbox.insert("end", sample)
            refresh()

        def from_file() -> None:
            value = filedialog.askopenfilename(filetypes=TSV)
            if not value:
                return
            wanted = {line.strip().split("\t")[0] for line in
                      Path(value).read_text(encoding="utf-8-sig", errors="replace").splitlines() if line.strip()}
            items = listbox.get(0, "end")
            for n, sample in enumerate(items):
                if sample in wanted:
                    listbox.selection_set(n)
            missing = wanted - set(items)
            refresh()
            if missing:
                messagebox.showwarning("Varan", f"{len(missing)} sample(s) of the file are not in the study:\n"
                                       + ", ".join(sorted(missing)[:20]))

        study.trace_add("write", load)
        listbox.bind("<<ListboxSelect>>", refresh)
        box = section(body, "Output", 2)
        self._output_fields(box, 0, operation)
        ttk.Checkbutton(box, text="Overwrite the last version", variable=overwrite).grid(
            row=3, column=1, sticky="w", pady=(6, 0))
        setattr(self, f"{operation}_vars", (study, listbox, overwrite))
        self._buttons(body, 3, run=lambda: self.start(self._study_form(operation)))

    def _build_settings(self) -> None:
        page = Scrollable(self.tabs)
        self.tabs.add(page, text="Settings")
        body = page.body
        body.columnconfigure(0, weight=1)
        s = self.settings
        self.set_vars = {k: tk.StringVar(value=getattr(s, k)) for k in
                         ("engine", "image", "vep_cache", "vep_version", "ref_fasta", "oncokb_token", "conf_ini",
                          "intro")}
        v = self.set_vars

        box = section(body, "Container engine", 0)
        ttk.Label(box, text="Engine").grid(row=0, column=0, sticky="w", padx=(0, 12), pady=3)
        engine = ttk.Combobox(box, textvariable=v["engine"], state="readonly", values=["docker", "podman"], width=12)
        engine.grid(row=0, column=1, sticky="w", pady=3)
        engine.bind("<<ComboboxSelected>>", lambda e: self.check_engine(start=True))
        ttk.Label(box, text="Varan image").grid(row=1, column=0, sticky="w", padx=(0, 12), pady=3)
        ttk.Entry(box, textvariable=v["image"], width=60).grid(row=1, column=1, sticky="ew", pady=3)
        ttk.Button(box, text="Download / update", command=self.pull).grid(row=1, column=2, padx=(6, 0))
        self.image_label = ttk.Label(box, text="", foreground="gray40")
        self.image_label.grid(row=2, column=1, sticky="w")

        box = section(body, "Reference data (needed for new studies)", 1)
        PathField(box, 0, "VEP cache folder", v["vep_cache"], kind="dir", app=self)
        ttk.Label(box, text="VEP cache version").grid(row=1, column=0, sticky="w", padx=(0, 12), pady=3)
        self.version_box = ttk.Combobox(box, textvariable=v["vep_version"], width=8)
        self.version_box.grid(row=1, column=1, sticky="w", pady=3)
        PathField(box, 2, "Reference FASTA (GRCh37)", v["ref_fasta"],
                  filetypes=[("FASTA", "*.fa *.fasta *.fa.gz *.fna"), ("All files", "*.*")], app=self)
        ttk.Label(box, text="Cache: the folder with homo_sapiens/<version>_GRCh37 inside.\n"
                  "FASTA: the genome, with its .fai index next to it.",
                  foreground="gray40").grid(row=3, column=1, sticky="w")
        v["vep_cache"].trace_add("write", lambda *_: self._vep_versions())
        self._vep_versions()

        box = section(body, "OncoKB", 2)
        ttk.Label(box, text="API token").grid(row=0, column=0, sticky="w", padx=(0, 12), pady=3)
        token = ttk.Entry(box, textvariable=v["oncokb_token"], show="•", width=50)
        token.grid(row=0, column=1, sticky="w", pady=3)
        link = ttk.Label(box, text="Get a token at oncokb.org", foreground="#1a6fb3", cursor="hand2")
        link.grid(row=1, column=1, sticky="w")
        link.bind("<Button-1>", lambda e: webbrowser.open("https://www.oncokb.org/account/register"))

        box = section(body, "Advanced", 3)
        PathField(box, 0, "Own conf.ini", v["conf_ini"], filetypes=[("conf.ini", "*.ini"), ("All files", "*.*")],
                  app=self, optional=True)
        ttk.Label(box, text="Its values replace Varan's defaults (filters, thresholds, cBioPortal headers...).\n"
                  "Paths are set by this window.", foreground="gray40").grid(row=1, column=1, sticky="w")
        tools = ttk.Frame(box)
        tools.grid(row=2, column=1, sticky="w", pady=(6, 0))
        ttk.Button(tools, text="Make a conf.ini to edit...", command=self.export_conf).pack(side="left", padx=(0, 8))
        ttk.Button(tools, text="Open the runs folder", command=lambda: open_path(app_dir() / "runs")).pack(side="left")
        ttk.Checkbutton(box, text="Show the intro video at start", variable=v["intro"], onvalue="yes",
                        offvalue="no", command=self._save_settings).grid(row=3, column=1, sticky="w", pady=(8, 0))
        ttk.Label(body, text=f"Varan app {__version__}", foreground="gray50").grid(row=4, column=0, sticky="e")

    def _build_log(self, parent) -> ttk.Frame:
        frame = ttk.Frame(parent, padding=(0, 8, 0, 8))
        status = ttk.Frame(frame)
        status.pack(fill="x")
        self.stage_label = ttk.Label(status, text="Ready.")
        self.stage_label.pack(side="left")
        self.progress = ttk.Progressbar(status, mode="indeterminate", length=180)
        self.stop_button = ttk.Button(status, text="Stop", command=self.stop, state="disabled")
        self.stop_button.pack(side="right")
        self.folder_button = ttk.Button(status, text="Open study folder", state="disabled",
                                        command=lambda: self.result and open_path(self.result))
        self.folder_button.pack(side="right", padx=6)
        self.report_button = ttk.Button(status, text="Open report", state="disabled", command=self.open_report)
        self.report_button.pack(side="right")
        self.log = ScrolledText(frame, height=8, wrap="word", font=("Consolas", 9) if sys.platform == "win32"
                                else ("Menlo", 11) if sys.platform == "darwin" else ("monospace", 9))
        self.log.pack(fill="both", expand=True, pady=(6, 0))
        self.log.tag_configure("error", foreground="#b3261e")
        self.log.tag_configure("warning", foreground="#a15c00")
        self.log.tag_configure("success", foreground="#1e7b34")
        self.log.tag_configure("app", foreground="#1a6fb3")
        self.log.configure(state="disabled")
        return frame

    # ------------------------------------------------------- form helpers

    def _input_mode(self) -> None:
        table = self.create_vars["input_mode"].get() == "table"
        for field in self.table_fields:
            field.enable(table)
        self.folder_field.enable(not table)

    def _pipeline(self) -> None:
        key = self._pipeline_key()
        if key in jobs.IMPLIED_TYPE:
            self.create_vars["sample_type"].set(jobs.IMPLIED_TYPE[key])
            self.type_box.configure(state="disabled")
        else:
            self.type_box.configure(state="readonly")

    def _pipeline_key(self) -> str:
        label = self.create_vars["pipeline"].get()
        return next((k for k, v in jobs.PIPELINES.items() if v == label), "")

    def _filter_cancers(self, box: ttk.Combobox) -> None:
        text = box.get().strip().lower()
        box.configure(values=[c for c in self.cancers if text in c.lower()] if text else self.cancers)

    def _vep_versions(self) -> None:
        folder = self.set_vars["vep_cache"].get().strip()
        versions = vep_cache_versions(folder) if folder else []
        self.version_box.configure(values=versions)
        if versions and self.set_vars["vep_version"].get() not in versions:
            self.set_vars["vep_version"].set(versions[0])

    def _save_settings(self) -> None:
        for key, var in self.set_vars.items():
            setattr(self.settings, key, var.get().strip())
        try:
            self.settings.save()
        except OSError as err:
            self.say(f"Settings not saved: {err}", "warning")

    def _create_form(self) -> jobs.CreateForm:
        v = self.create_vars
        cancer = v["cancer"].get().strip()
        return jobs.CreateForm(
            input_mode=v["input_mode"].get(), sample_tsv=v["sample_tsv"].get(), patient_tsv=v["patient_tsv"].get(),
            fusion_tsv=v["fusion_tsv"].get(), input_folder=v["input_folder"].get(),
            pipeline=self._pipeline_key(), sample_type=v["sample_type"].get(),
            cancer=cancer.split(" - ")[0].strip(), project_id=v["project_id"].get(),
            project_name=v["project_name"].get(), description=v["description"].get(),
            output_dir=self.create_out_dir.get(), study_name=self.create_name.get(),
            oncokb=v["oncokb"].get(), multiple=v["multiple"].get(), sigma=v["sigma"].get(),
            filters="".join(k for k, var in self.filter_vars.items() if var.get()),
            existing=v["existing"].get())

    def _study_form(self, operation: str) -> jobs.StudyForm:
        if operation == "update":
            return jobs.StudyForm("update", self.update_study.get(), self.update_new.get(), [],
                                  self.update_out_dir.get(), self.update_name.get(), self.update_overwrite.get())
        study, listbox, overwrite = getattr(self, f"{operation}_vars")
        samples = [listbox.get(i) for i in listbox.curselection()]
        return jobs.StudyForm(operation, study.get(), "", samples, getattr(self, f"{operation}_out_dir").get(),
                              getattr(self, f"{operation}_name").get(), overwrite.get())

    def export_conf(self) -> None:
        dest = filedialog.asksaveasfilename(defaultextension=".ini", initialfile="conf.ini",
                                            filetypes=[("conf.ini", "*.ini")])
        if not dest:
            return
        Path(dest).write_text(resource("default_conf.ini").read_text(encoding="utf-8"), encoding="utf-8")
        self.set_vars["conf_ini"].set(str(Path(dest)))
        self._save_settings()
        open_path(dest)

    # ------------------------------------------------------------ the log

    def say(self, text: str, tag: str = "app") -> None:
        self.events.put(("log", text, tag))

    def _append(self, text: str, tag: str | None) -> None:
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n", tag or ())
        self.log.see("end")
        self.log.configure(state="disabled")

    def _line(self, line: str) -> None:
        line = ANSI.sub("", line)
        if self.run is not None:
            line = self.run.host_text(line)
        tag = None
        if "| CRITICAL" in line or "| ERROR" in line or line.startswith("docker:"):
            tag = "error"
        elif "| WARNING" in line:
            tag = "warning"
        elif "| SUCCESS" in line:
            tag = "success"
        for marker, stage in STAGES:
            if marker in line:
                if marker.startswith("Starting vcf2maf"):
                    sample = line.split(marker, 1)[1].split("...")[0].strip()
                    stage = f"{stage}: {sample}"
                self.events.put(("stage", stage))
        if self.run is not None:
            try:
                with (self.run.staging / "varan_gui.log").open("a", encoding="utf-8") as handle:
                    handle.write(line + "\n")
            except OSError:
                pass
        self.events.put(("log", line, tag))

    def _poll(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                kind = event[0]
                if kind == "log":
                    self._append(event[1], event[2])
                elif kind == "stage":
                    self.stage_label.configure(text=event[1])
                elif kind == "engine":
                    self._show_engine(*event[1:])
                elif kind == "image":
                    self.image_label.configure(text=event[1])
                elif kind == "done":
                    self._finished(*event[1:])
                elif kind == "call":
                    event[1]()
        except queue.Empty:
            pass
        self.after(100, self._poll)

    # ------------------------------------------------------------- engine

    def check_engine(self, start: bool = False) -> None:
        self._save_settings() if hasattr(self, "set_vars") else None
        self.engine = Engine(self.settings.engine)
        engine = self.engine

        def work() -> None:
            name = "Docker" if engine.kind == "docker" else "Podman"
            if not engine.installed:
                self.events.put(("engine", f"{name} is not installed", "install"))
                return
            if not engine.running():
                if not start or not engine.start():
                    self.events.put(("engine", f"{name} is not running - start it, then Retry", "retry"))
                    return
                self.events.put(("engine", f"Starting {name}...", None))
                if not engine.wait():
                    self.events.put(("engine", f"{name} did not start - start it, then Retry", "retry"))
                    return
            self.events.put(("engine", f"{name} ready", "ready"))
            self._image_status()

        threading.Thread(target=work, daemon=True).start()

    def _show_engine(self, text: str, state: str | None) -> None:
        self.engine_ready = state == "ready"
        self.engine_label.configure(text=("● " if self.engine_ready else "○ ") + text,
                                    foreground="#1e7b34" if self.engine_ready else "#a15c00")
        if state == "install":
            self.engine_button.configure(text="Get it", command=lambda: webbrowser.open(DOWNLOAD_PAGE[self.engine.kind]))
            self.engine_button.pack(side="left")
        elif state == "retry":
            self.engine_button.configure(text="Retry", command=lambda: self.check_engine(start=True))
            self.engine_button.pack(side="left")
        else:
            self.engine_button.pack_forget()

    def _image_status(self) -> bool:
        image = self.settings.image
        found = bool(self.engine.image_id(image))
        self.events.put(("image", "Downloaded." if found else
                         "Not downloaded yet: it is fetched on the first run (about 3 GB)."))
        return found

    def pull(self) -> None:
        if self.busy:
            return
        self._save_settings()
        if not self.engine_ready:
            messagebox.showwarning("Varan", "Docker is not running yet.")
            return
        self._begin(f"Downloading {self.settings.image}...")
        threading.Thread(target=lambda: self._pull_then(None), daemon=True).start()

    def _pull_then(self, after) -> None:
        self.say(f"Downloading the Varan image {self.settings.image} - only the first time, "
                 "it can take a while.")
        code = self.engine.stream(["pull", *self.engine.platform_args(), self.settings.image],
                                  lambda line: self.events.put(("log", ANSI.sub("", line), None)),
                                  on_start=self._started)
        self._image_status()
        if code != 0:
            self.events.put(("done", code, None, "The image could not be downloaded."))
        elif after is None:
            self.events.put(("done", 0, None, "Image downloaded."))
        else:
            after()

    # ---------------------------------------------------------------- run

    def start(self, form, check: bool = False) -> None:
        if self.busy:
            return
        self._save_settings()
        if not self.engine_ready:
            messagebox.showwarning("Varan", "Docker is not running yet: wait for the green light at the top right.")
            return
        try:
            run = jobs.prepare(form, self.settings, self.engine.platform_args())
        except jobs.FormError as err:
            messagebox.showerror("Varan", str(err))
            return
        except OSError as err:
            messagebox.showerror("Varan", f"Could not prepare the run: {err}")
            return
        if check:
            run.args.append("-D")
        if run.warnings and not messagebox.askokcancel("Varan", "\n\n".join(run.warnings) + "\n\nRun anyway?"):
            return
        self.run, self.result = run, None
        self.report_button.configure(state="disabled")
        self.folder_button.configure(state="disabled")
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")
        self._begin("Checking the inputs..." if check else "Starting...")
        self.say(jobs.command_text(run.args))

        def work() -> None:
            if not self.engine.image_id(self.settings.image):
                self._pull_then(lambda: self._execute(run, check))
            else:
                self._execute(run, check)

        threading.Thread(target=work, daemon=True).start()

    def _execute(self, run: jobs.Run, check: bool) -> None:
        try:
            code = self.engine.stream(run.args, self._line, on_start=self._started)
        except OSError as err:
            self.say(str(err), "error")
            code = -1
        result = None if check else jobs.find_result(run)
        if result is not None:
            jobs.write_record(run, result, self.engine.image_id(self.settings.image), __version__)
        if check:
            message = "The inputs look fine." if code == 0 else "Some inputs need fixing: see the log."
        elif self.stopping:
            message = "Stopped." + (' To go on later: "resume the last version".' if result else "")
        elif code == 0:
            message = "Done."
        else:
            message = "Varan stopped with an error: see the log" + (" and the report." if result else ".")
        self.events.put(("done", code, result, message))

    def _started(self, proc: subprocess.Popen) -> None:
        self.proc = proc

    def _begin(self, text: str) -> None:
        self.busy, self.stopping = True, False
        self.stage_label.configure(text=text)
        self.progress.pack(side="left", padx=10, after=self.stage_label)
        self.progress.start(12)
        self.stop_button.configure(state="normal")
        for button in self.run_buttons:
            button.configure(state="disabled")

    def _finished(self, code: int, result: Path | None, message: str) -> None:
        self.busy, self.proc = False, None
        self.progress.stop()
        self.progress.pack_forget()
        self.stop_button.configure(state="disabled")
        for button in self.run_buttons:
            button.configure(state="normal")
        self.stage_label.configure(text=message)
        self.result = result
        if result is not None:
            self.folder_button.configure(state="normal")
            if (result / "report_VARAN.html").exists():
                self.report_button.configure(state="normal")
            self._append(f"Study folder: {result}", "success" if code == 0 else "warning")
        if code == 0 and result is not None and (result / "report_VARAN.html").exists():
            self.open_report()
        elif code != 0 and not self.stopping and self.run is not None:
            self._append(f"Full log: {self.run.staging / 'varan_gui.log'}", "app")

    def open_report(self) -> None:
        if self.result is not None and (self.result / "report_VARAN.html").exists():
            open_path(self.result / "report_VARAN.html")

    def stop(self) -> None:
        if not self.busy:
            return
        self.stopping = True
        self.stage_label.configure(text="Stopping...")
        name = self.run.name if self.run is not None else None

        def work() -> None:
            if name:
                self.engine.kill(name)
            if self.proc is not None and self.proc.poll() is None:
                self.proc.terminate()

        threading.Thread(target=work, daemon=True).start()

    def close(self) -> None:
        if self.busy and not messagebox.askyesno("Varan", "Varan is still running. Stop it and quit?"):
            return
        if self.busy:
            self.stop()
        self._save_settings()
        self.destroy()


def main() -> None:
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except (AttributeError, OSError):
            pass
    App().mainloop()


if __name__ == "__main__":
    main()
