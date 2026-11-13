"""Turn assets/intro.mp4 into the frames the app plays at start.

tkinter can't play videos, but decodes a GIF image in about a millisecond:
each frame becomes a GIF in assets/intro_frames.zip (with fps.txt), which
build.py bundles in the app. Not in git - this makes it again:

    pip install pillow imageio imageio-ffmpeg
    python -m varan_gui.make_intro
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

ASSETS = Path(__file__).resolve().parent / "assets"
VIDEO = ASSETS / "intro.mp4"
FRAMES = ASSETS / "intro_frames.zip"
WIDTH = 640
STEP = 2  # every other frame: 24 fps -> 12


def make(video: Path = VIDEO, dest: Path = FRAMES) -> Path:
    import imageio.v3 as iio
    from PIL import Image

    fps = iio.immeta(video).get("fps", 24) / STEP
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_STORED) as archive:
        archive.writestr("fps.txt", f"{fps:g}")
        for n, frame in enumerate(iio.imiter(video)):
            if n % STEP:
                continue
            image = Image.fromarray(frame)
            image = image.resize((WIDTH, round(image.height * WIDTH / image.width)), Image.LANCZOS)
            buf = io.BytesIO()
            image.quantize(colors=256, method=Image.Quantize.MEDIANCUT,
                           dither=Image.Dither.FLOYDSTEINBERG).save(buf, "GIF", optimize=True)
            archive.writestr(f"{n // STEP:04d}.gif", buf.getvalue())
    return dest


def up_to_date() -> bool:
    return FRAMES.exists() and (not VIDEO.exists() or FRAMES.stat().st_mtime >= VIDEO.stat().st_mtime)


if __name__ == "__main__":
    print(make())
