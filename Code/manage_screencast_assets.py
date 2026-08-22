"""
manage_screencast_assets.py
-----------------------------------------------------------------------------
Reusable utility for populating Docs/screencasts/ (git-ignored -- see
.gitignore -- so raw video/image assets never get committed; only the
YouTube links embedded in User_Guide.md do).

Two modes:

1. Build/extend the intro clip from one or more still images, each held for
   a fixed duration (default 1s), concatenated into Docs/screencasts/00_intro.mp4:

       python manage_screencast_assets.py --intro "C:\\path\\a.jpg" "C:\\path\\b.png"

   Re-running with --intro replaces 00_intro.mp4 (rebuilds from the images
   given on that invocation, in the order given).

2. Add recorded screencast clips to the folder with the next available
   two-digit sequence prefix (01_, 02_, ...), so they sort after the intro
   and before Youtube upload:

       python manage_screencast_assets.py --add "C:\\path\\step1.mp4" "C:\\path\\step2.mp4"

3. Prepend 00_intro.mp4 onto one or more recorded clips (re-encoding the
   intro to match each recording's resolution/fps and adding a silent audio
   track so the streams concatenate cleanly), writing each result into
   Docs/screencasts/ with the next sequence prefix. Originals are left
   untouched -- only new, merged files are created:

       python manage_screencast_assets.py --merge-intro "C:\\path\\step1.mp4" "C:\\path\\step2.mp4"

Requires: pillow, imageio, imageio-ffmpeg (installed into the system Python
via `python -m pip install pillow imageio imageio-ffmpeg`).
-----------------------------------------------------------------------------
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image
import imageio.v2 as imageio
import imageio_ffmpeg
import numpy as np

SCREENCASTS_DIR = Path(__file__).resolve().parent.parent / "Docs" / "screencasts"
INTRO_NAME = "00_intro.mp4"
FPS = 25
FFMPEG = imageio_ffmpeg.get_ffmpeg_exe()


def build_intro(image_paths, seconds_per_image=1.0, canvas=(720, 720), bg=(255, 255, 255)):
    SCREENCASTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = SCREENCASTS_DIR / INTRO_NAME

    frames_per_image = max(1, round(seconds_per_image * FPS))
    writer = imageio.get_writer(str(out_path), fps=FPS, codec="libx264",
                                 quality=8, macro_block_size=1)
    try:
        for p in image_paths:
            src = Path(p)
            if not src.exists():
                print(f"  SKIP (not found): {src}")
                continue
            img = Image.open(src).convert("RGB")
            # Letterbox onto a fixed canvas so mismatched image sizes don't
            # break the video stream (all frames must share one resolution).
            frame = Image.new("RGB", canvas, bg)
            scale = min(canvas[0] / img.width, canvas[1] / img.height)
            new_size = (max(1, round(img.width * scale)), max(1, round(img.height * scale)))
            resized = img.resize(new_size, Image.LANCZOS)
            offset = ((canvas[0] - new_size[0]) // 2, (canvas[1] - new_size[1]) // 2)
            frame.paste(resized, offset)
            arr = np.array(frame)
            for _ in range(frames_per_image):
                writer.append_data(arr)
            print(f"  Added: {src.name}  ({seconds_per_image}s)")
    finally:
        writer.close()

    print(f"\nIntro clip written: {out_path}  "
          f"({len(image_paths)} image(s) x {seconds_per_image}s)")


def _next_sequence_number():
    import re
    existing = [int(m.group(1)) for f in SCREENCASTS_DIR.glob("*")
                if (m := re.match(r"^(\d{2})_", f.name))]
    return (max(existing) + 1) if existing else 1  # 00 is reserved for the intro


def _probe(path):
    """Return (width, height, fps, duration_seconds) parsed from `ffmpeg -i`."""
    import re
    proc = subprocess.run([FFMPEG, "-i", str(path)], capture_output=True, text=True)
    text = proc.stderr
    res_m = re.search(r"Video:.*?(\d{2,5})x(\d{2,5})", text)
    width, height = int(res_m.group(1)), int(res_m.group(2))
    fps_m = re.search(r"([\d.]+)\s*fps", text)
    fps = float(fps_m.group(1)) if fps_m else 30.0
    dur_m = re.search(r"Duration:\s*(\d+):(\d+):([\d.]+)", text)
    duration = (int(dur_m.group(1)) * 3600 + int(dur_m.group(2)) * 60 + float(dur_m.group(3))
                if dur_m else 2.0)
    return width, height, fps, duration


def add_clips(video_paths):
    SCREENCASTS_DIR.mkdir(parents=True, exist_ok=True)
    next_n = _next_sequence_number()

    for p in video_paths:
        src = Path(p)
        if not src.exists():
            print(f"  SKIP (not found): {src}")
            continue
        dest = SCREENCASTS_DIR / f"{next_n:02d}_{src.name}"
        shutil.copy2(src, dest)
        print(f"  Added: {dest.name}")
        next_n += 1


def merge_intro(video_paths):
    """Prepend 00_intro.mp4 onto each recording, writing a new sequenced file.

    The intro is re-scaled/padded to each recording's own resolution and fps,
    and given a silent audio track matching its duration, so the ffmpeg concat
    filter can join heterogeneous inputs (the intro is a plain 720x720/25fps
    synthetic clip with no audio; recordings are screen captures with audio at
    whatever resolution/fps OBS or similar produced). Originals are untouched.
    """
    SCREENCASTS_DIR.mkdir(parents=True, exist_ok=True)
    intro_path = SCREENCASTS_DIR / INTRO_NAME
    if not intro_path.exists():
        print(f"  ERROR: {INTRO_NAME} not found in {SCREENCASTS_DIR} -- build it first with --intro.")
        return

    _, _, _, intro_duration = _probe(intro_path)
    next_n = _next_sequence_number()

    for p in video_paths:
        src = Path(p)
        if not src.exists():
            print(f"  SKIP (not found): {src}")
            continue

        width, height, fps, _ = _probe(src)
        dest = SCREENCASTS_DIR / f"{next_n:02d}_{src.name}"

        filter_complex = (
            f"[0:v]scale={width}:{height}:force_original_aspect_ratio=decrease,"
            f"pad={width}:{height}:(ow-iw)/2:(oh-ih)/2,setsar=1,fps={fps},format=yuv420p[v0];"
            f"[1:a]aformat=sample_rates=48000:channel_layouts=stereo[a0];"
            f"[2:v]fps={fps},format=yuv420p[v1];"
            f"[2:a]aformat=sample_rates=48000:channel_layouts=stereo[a1];"
            f"[v0][a0][v1][a1]concat=n=2:v=1:a=1[outv][outa]"
        )
        cmd = [
            FFMPEG, "-y",
            "-i", str(intro_path),
            "-f", "lavfi", "-t", str(intro_duration), "-i", "anullsrc=r=48000:cl=stereo",
            "-i", str(src),
            "-filter_complex", filter_complex,
            "-map", "[outv]", "-map", "[outa]",
            "-c:v", "libx264", "-c:a", "aac", "-movflags", "+faststart",
            str(dest),
        ]
        print(f"  Merging intro + {src.name} -> {dest.name} ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  FFMPEG ERROR for {src.name}:\n{result.stderr[-2000:]}")
            continue
        print(f"  Done: {dest.name}")
        next_n += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--intro", nargs="+", metavar="IMAGE",
                     help="Build Docs/screencasts/00_intro.mp4 from these images, in order.")
    ap.add_argument("--seconds", type=float, default=1.0,
                     help="Seconds each image is held in the intro clip (default: 1.0).")
    ap.add_argument("--add", nargs="+", metavar="VIDEO",
                     help="Copy these recorded clips into Docs/screencasts/ with the next "
                          "sequence prefix.")
    ap.add_argument("--merge-intro", nargs="+", metavar="VIDEO",
                     help="Prepend 00_intro.mp4 onto each of these recordings and write the "
                          "result into Docs/screencasts/ with the next sequence prefix.")
    args = ap.parse_args()

    if not args.intro and not args.add and not args.merge_intro:
        ap.print_help()
        sys.exit(1)

    if args.intro:
        build_intro(args.intro, seconds_per_image=args.seconds)
    if args.add:
        add_clips(args.add)
    if args.merge_intro:
        merge_intro(args.merge_intro)


if __name__ == "__main__":
    main()
