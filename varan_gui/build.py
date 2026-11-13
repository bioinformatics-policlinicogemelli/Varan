"""Build the double-clickable app with PyInstaller.

    pip install pyinstaller pillow imageio imageio-ffmpeg
    python -m varan_gui.build

Windows and Linux get dist/Varan(.exe), one file; macOS gets dist/Varan.app.
Run it on each system the app is for: PyInstaller doesn't cross-compile.
"""

from __future__ import annotations

import sys
from pathlib import Path

import PyInstaller.__main__

from varan_gui import make_intro

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    sep = ";" if sys.platform == "win32" else ":"
    gui = ROOT / "varan_gui"
    assets = gui / "assets"
    icon = assets / ("varan.ico" if sys.platform == "win32" else "varan.png")
    if not make_intro.up_to_date():
        print("Making the intro frames from", make_intro.VIDEO)
        make_intro.make()
    args = [
        str(gui / "__main__.py"),
        "--name", "Varan",
        "--windowed",
        "--noconfirm",
        "--clean",
        "--icon", str(icon),
        "--paths", str(ROOT),
        *(arg for name in ("varan.png", "varan_header.png", "intro_frames.zip")
          for arg in ("--add-data", f"{assets / name}{sep}varan_gui/assets")),
        "--add-data", f"{gui / 'default_conf.ini'}{sep}varan_gui",
        "--add-data", f"{ROOT / 'cancer_list.txt'}{sep}.",
        "--distpath", str(ROOT / "dist"),
        "--workpath", str(ROOT / "build"),
        "--specpath", str(ROOT / "build"),
    ]
    if sys.platform == "darwin":
        args += ["--osx-bundle-identifier", "org.policlinicogemelli.varan"]
    else:
        args.append("--onefile")
    PyInstaller.__main__.run(args)


if __name__ == "__main__":
    main()
