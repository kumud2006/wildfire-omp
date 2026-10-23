#!/usr/bin/env python3
"""
viz.py -- convert the PPM images written by ./wildfire into PNG (needs Pillow).

  python3 viz.py fire.ppm                      -> fire.png
  python3 viz.py fire_t00*.ppm --montage m.png -> every file as PNG + one side-by-side strip

Colours: green = unburnt (brighter = more fuel), orange = burning, grey = burned.
"""
import argparse
import os
from PIL import Image, ImageDraw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--montage", help="also write all frames side by side to this PNG")
    ap.add_argument("--max", type=int, default=600, help="max width/height of each frame in pixels")
    a = ap.parse_args()

    frames = []
    for f in a.files:
        im = Image.open(f).convert("RGB")
        if max(im.size) > a.max:                       # downscale big grids for viewing
            im.thumbnail((a.max, a.max), Image.NEAREST)
        out = os.path.splitext(f)[0] + ".png"
        im.save(out)
        print("wrote", out)
        frames.append((f, im))

    if a.montage and len(frames) > 1:
        w = sum(im.width for _, im in frames) + 8 * (len(frames) - 1)
        h = max(im.height for _, im in frames) + 18
        sheet = Image.new("RGB", (w, h), (255, 255, 255))
        d = ImageDraw.Draw(sheet)
        x = 0
        for f, im in frames:
            sheet.paste(im, (x, 18))
            d.text((x + 2, 3), os.path.basename(f), fill=(0, 0, 0))
            x += im.width + 8
        sheet.save(a.montage)
        print("wrote", a.montage)


if __name__ == "__main__":
    main()
