"""What the user filled in, turned into one `docker run` of Varan.

`prepare()` checks the form, writes the run's staging folder (conf.ini and
the rewritten tables) and returns the engine arguments; `find_result()`
finds the study folder the run made.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from varan_gui import staging
from varan_gui.settings import Settings, app_dir, cancer_types, resource, vep_cache_root, vep_cache_versions

PIPELINES = {
    "": "Generic (sample.tsv with VCF files)",
    "illumina_solid": "Illumina TSO500 - tissue",
    "illumina_liquid": "Illumina TSO500 - plasma",
    "guardant": "Guardant360",
}
# Pipelines that fix the sample type (varan.py's implied_sample_type()).
IMPLIED_TYPE = {"illumina_solid": "Solid", "illumina_liquid": "Liquid", "guardant": "Liquid"}

FILTERS = {
    "d": "Drop calls without ALT or not PASS (raw VCF)",
    "p": "Keep only FILTER = PASS (MAF)",
    "v": "Tumour VAF between t_VAF_min and t_VAF_max",
    "n": "Novel variants need t_VAF_min_novel (with v)",
    "a": "Population allele frequency (AF)",
    "c": "Drop benign ClinVar significance",
    "q": "Keep the listed consequences",
    "y": "PolyPhen",
    "i": "Drop LOW / MODIFIER impact",
    "o": "OncoKB oncogenic only (needs OncoKB)",
}

KEEP_RUNS = 30


@dataclass
class CreateForm:
    input_mode: str = "table"  # table or folder
    sample_tsv: str = ""
    patient_tsv: str = ""
    fusion_tsv: str = ""
    input_folder: str = ""
    pipeline: str = ""
    sample_type: str = "Solid"
    cancer: str = ""
    project_id: str = ""
    project_name: str = ""
    description: str = ""
    output_dir: str = ""
    study_name: str = ""
    oncokb: bool = False
    multiple: bool = False
    sigma: bool = False
    filters: str = ""
    existing: str = "new"  # new, overwrite or resume


@dataclass
class StudyForm:
    operation: str = "update"  # update, extract or remove
    study: str = ""
    new_study: str = ""
    samples: list[str] = field(default_factory=list)
    output_dir: str = ""
    study_name: str = ""
    overwrite: bool = False


@dataclass
class Run:
    staging: Path
    name: str
    args: list[str]
    output_dir: Path
    output_base: str
    started: float
    warnings: list[str]
    paths: list[tuple[str, Path]] = field(default_factory=list)

    def host_text(self, text: str) -> str:
        """`text` with container paths turned back into the user's paths."""
        def host(match: re.Match) -> str:
            for prefix, folder in self.paths:
                if match.group(0) == prefix or match.group(0).startswith(prefix + "/"):
                    rest = match.group(0)[len(prefix):].strip("/")
                    return str(folder.joinpath(*rest.split("/"))) if rest else str(folder)
            return match.group(0)
        return CONTAINER_PATH.sub(host, text)


CONTAINER_PATH = re.compile(r"/varan/(?:run|data/\d+)(?:/[^\s'\"\],;)]*)?")


def timezone_env() -> list[str]:
    """The computer's UTC offset as a POSIX TZ, so that the times in the
    log and the report are local ones (the image has no time zones)."""
    offset = datetime.now().astimezone().utcoffset()
    minutes = int(offset.total_seconds() // 60) if offset else 0
    if minutes == 0:
        return ["-e", "TZ=UTC0"]
    sign, minutes = ("+", minutes) if minutes > 0 else ("-", -minutes)
    hh, mm = divmod(minutes, 60)
    posix = "-" if sign == "+" else "+"
    return ["-e", f"TZ=<{sign}{hh:02d}{mm:02d}>{posix}{hh:02d}:{mm:02d}"]


class FormError(ValueError):
    """A problem the user has to fix before the run starts."""


def _existing(value: str, what: str, folder: bool = False) -> Path:
    if not value.strip():
        raise FormError(f"Choose {what}.")
    path = Path(value.strip())
    if not (path.is_dir() if folder else path.is_file()):
        raise FormError(f"{what[0].upper() + what[1:]} not found: {path}")
    return path


def _study_name(value: str) -> str:
    name = value.strip()
    if not name:
        raise FormError("Give the study folder a name.")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise FormError(f'Study folder name "{name}": use letters, digits, "-", "_" or "." (no spaces).')
    return re.sub(r"_v\d+$", "", name)


def new_staging() -> Path:
    runs = app_dir() / "runs"
    runs.mkdir(exist_ok=True)
    old = sorted(p for p in runs.iterdir() if p.is_dir())
    for stale in old[:-KEEP_RUNS]:
        shutil.rmtree(stale, ignore_errors=True)
    folder = runs / datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    folder.mkdir()
    return folder


def _base_conf(settings: Settings) -> tuple:
    """Default conf.ini + the user's own + Settings; and the folder its
    relative paths are relative to."""
    conf = staging.read_conf(resource("default_conf.ini"))
    base = Path.cwd()
    if settings.conf_ini.strip():
        user = Path(settings.conf_ini.strip())
        if not user.is_file():
            raise FormError(f"The conf.ini chosen in Settings is missing: {user}")
        staging.merge_conf(conf, staging.read_conf(user))
        base = user.parent
    if settings.oncokb_token.strip():
        conf.set("OncoKB", "ONCOKB", settings.oncokb_token.strip())
    conf.set("Paths", "CACHE", settings.vep_version.strip() or conf.get("Paths", "CACHE", fallback="116") or "116")
    for key, value in staging.IMAGE_PATHS.items():
        conf.set("Paths", key, value)
    return conf, base


def _reference_data(conf, settings: Settings, mounts: staging.Mounts, warnings: list[str]) -> tuple:
    """VEP cache, FASTA and ClinVar for a create run."""
    cache = vep_cache_root(_existing(settings.vep_cache, "the VEP cache folder (Settings)", folder=True))
    if cache is None:
        raise FormError("The VEP cache folder (Settings) has no homo_sapiens folder inside.")
    version = settings.vep_version.strip()
    if version and version not in vep_cache_versions(cache):
        raise FormError(f"The VEP cache has no homo_sapiens/{version}_GRCh37 folder: "
                        "choose the version found in it (Settings).")
    fasta = _existing(settings.ref_fasta, "the reference FASTA (Settings)")
    if not Path(f"{fasta}.fai").exists():
        warnings.append(f"No index {fasta.name}.fai next to the FASTA: vcf2maf may fail.")
    mounts.add(cache)
    mounts.add(fasta.parent)
    clinv = conf.get("Paths", "CLINV", fallback="")
    clinv_file = staging.as_host_path(clinv.split(",")[0], Path.cwd()) if clinv.strip() else None
    if clinv_file is not None:
        mounts.add(clinv_file.parent)
    return cache, fasta, clinv_file


def prepare(form: CreateForm | StudyForm, settings: Settings, platform_args: list[str]) -> Run:
    warnings: list[str] = []
    mounts = staging.Mounts()
    conf, conf_base = _base_conf(settings)
    output_dir = _existing(form.output_dir, "the output folder", folder=True)
    name = _study_name(form.study_name)
    mounts.add(output_dir)

    if isinstance(form, CreateForm):
        if not form.cancer.strip():
            raise FormError("Choose the cancer type.")
        known = {cid.lower() for cid, _ in cancer_types()}
        if known and form.cancer.strip().lower() not in known:
            raise FormError(f'"{form.cancer.strip()}" is not a cancer type of cBioPortal: choose one from the list.')
        if form.oncokb and not conf.get("OncoKB", "ONCOKB", fallback="").strip():
            raise FormError("OncoKB annotation needs an OncoKB token (Settings).")
        if "o" in form.filters and not form.oncokb:
            raise FormError('The OncoKB filter needs "OncoKB annotation" on.')
        if "n" in form.filters and "v" not in form.filters:
            raise FormError("The novel-variant filter needs the VAF filter on.")
        cache, fasta, clinv_file = _reference_data(conf, settings, mounts, warnings)
        conf.set("Vendor", "PIPELINE", form.pipeline)
        conf.set("Sample_Type", "TYPE", "" if form.pipeline in IMPLIED_TYPE else form.sample_type)
        for key, value in (("PROJECT_ID", form.project_id), ("PROJECT_NAME", form.project_name),
                           ("DESCRIPTION", form.description)):
            if value.strip():
                conf.set("Project", key, value.strip())
        tables = []
        if form.input_mode == "folder":
            folder = _existing(form.input_folder, "the input folder", folder=True)
            mounts.add(folder.parent)
        else:
            sample = _existing(form.sample_tsv, "the sample table (sample.tsv)")
            missing = staging.missing_path_cells(sample)
            if missing:
                shown = "\n".join(missing[:5]) + ("\n..." if len(missing) > 5 else "")
                raise FormError(f"{len(missing)} file(s) listed in {sample.name} not found:\n{shown}")
            tables.append(("sample.tsv", sample))
            for label, value in (("patient.tsv", form.patient_tsv), ("fusion.tsv", form.fusion_tsv)):
                tables.append((label, _existing(value, f"the {label} file") if value.strip() else None))
            for _, table in tables:
                if table is not None:
                    for host in staging.path_cells(table):
                        mounts.add(host if host.is_dir() else host.parent)
    else:
        # Parents, so that the container sees the folders under their own
        # names (name_v1...).
        study = _existing(form.study, "the study folder", folder=True)
        mounts.add(study.parent)
        if form.operation == "update":
            new_study = _existing(form.new_study, "the study to add", folder=True)
            mounts.add(new_study.parent)
        elif not form.samples:
            raise FormError("Choose at least one sample.")
        cache = fasta = clinv_file = None

    extra = staging.conf_host_paths(conf, conf_base)
    for _, _, host in extra:
        mounts.add(host if host.is_dir() else host.parent)

    # Every folder is known: now the container paths.
    run_dir = new_staging()
    if cache is not None:
        conf.set("Paths", "VEP_DATA", mounts.container_path(cache))
        conf.set("Paths", "REF_FASTA", mounts.container_path(fasta))
        if clinv_file is not None:
            rest = conf.get("Paths", "CLINV").split(",", 1)
            conf.set("Paths", "CLINV", ",".join([mounts.container_path(clinv_file), *rest[1:]]))
        else:
            conf.set("Paths", "CLINV", staging.IMAGE_CLINVAR)
    for section, key, host in extra:
        conf.set(section, key, mounts.container_path(host))
    staging.write_conf(conf, run_dir / "conf.ini")

    out = f"{mounts.container_path(output_dir)}/{name}"
    run = staging.RUN_MOUNT
    if isinstance(form, CreateForm):
        if form.input_mode == "folder":
            inputs = [mounts.container_path(form.input_folder)]
        else:
            inputs = []
            for label, table in tables:
                if table is None:
                    inputs.append("")
                    continue
                staging.rewrite_table(table, mounts, run_dir / label)
                inputs.append(f"{run}/{label}")
            while inputs and inputs[-1] == "":
                inputs.pop()
        varan = ["-i", *inputs, "-o", out, "-c", form.cancer.strip()]
        varan += ["-k"] if form.oncokb else []
        varan += ["-m"] if form.multiple else []
        varan += ["-g"] if form.sigma else []
        varan += ["-f", form.filters] if form.filters else []
        varan += {"overwrite": ["-w"], "resume": ["-R"]}.get(form.existing, [])
    else:
        varan = ["-p", mounts.container_path(form.study), "-o", out]
        if form.operation == "update":
            varan = ["-u", *varan, "-n", mounts.container_path(form.new_study)]
        else:
            (run_dir / "samples.txt").write_text("\n".join(form.samples) + "\n", encoding="utf-8", newline="\n")
            varan = ["-e" if form.operation == "extract" else "-r", *varan, "-s", f"{run}/samples.txt"]
        varan += ["-w"] if form.overwrite else []

    container = f"varan-gui-{run_dir.name}"
    args = ["run", "--rm", "--name", container, *platform_args, *timezone_env(), *mounts.args(run_dir),
            settings.image.strip(), "-C", f"{run}/conf.ini", *varan]
    (run_dir / "command.txt").write_text(command_text(args), encoding="utf-8")
    paths = [(run, run_dir)] + [(f"{staging.DATA_ROOT}/{n}", f) for n, f in enumerate(mounts.folders)]
    return Run(run_dir, container, args, output_dir, name, time.time(), warnings, paths)


def command_text(args: list[str]) -> str:
    """The command as it could be typed, for the log (the OncoKB token is
    in conf.ini, never on the command line)."""
    return subprocess.list2cmdline(["docker", *args])


def find_result(run: Run) -> Path | None:
    """The newest <name>_v<n> study folder this run wrote."""
    best: tuple[int, Path] | None = None
    pattern = re.compile(rf"^{re.escape(run.output_base)}_v(\d+)$")
    for child in run.output_dir.iterdir():
        match = pattern.match(child.name)
        if not match or not child.is_dir():
            continue
        if child.stat().st_mtime < run.started - 5:
            report = child / "report_VARAN.html"
            if not report.exists() or report.stat().st_mtime < run.started - 5:
                continue
        version = int(match.group(1))
        if best is None or version > best[0]:
            best = (version, child)
    return best[1] if best else None


def write_record(run: Run, result: Path, image_id: str, version: str) -> None:
    """varan_app.txt in the study: how it was made, and where the
    /varan/... paths of its report and log are on this computer."""
    lines = [
        f"Made with the Varan app {version}, {datetime.fromtimestamp(run.started):%Y-%m-%d %H:%M}.",
        f"Image: {run.args[run.args.index('-C') - 1]} ({image_id or 'id unknown'})",
        "",
        "Command:",
        command_text(run.args),
        "",
        "Paths in the report and the log:",
        *(f"  {prefix} = {folder}" for prefix, folder in run.paths if prefix != staging.RUN_MOUNT),
        f"  {staging.RUN_MOUNT} = the app's run folder (conf.ini and inputs as passed to Varan): {run.staging}",
    ]
    try:
        (result / "varan_app.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    except OSError:
        pass


def study_samples(study: Path) -> list[str]:
    """SAMPLE_ID column of a study's data_clinical_sample.txt."""
    path = study / "data_clinical_sample.txt"
    if not path.is_file():
        return []
    rows = [line.split("\t") for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip() and not line.startswith("#")]
    if not rows or "SAMPLE_ID" not in rows[0]:
        return []
    col = rows[0].index("SAMPLE_ID")
    return [r[col].strip() for r in rows[1:] if len(r) > col and r[col].strip()]
