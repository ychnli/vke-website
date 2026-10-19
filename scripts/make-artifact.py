"""Build a copy of the site that can be published as a claude.ai Artifact.

The artifact host differs from a normal static host in four ways:
  * it adds its own <html>/<head>, so the page is written without them (title first);
  * paths must be relative, and no published path may start with "_" (Astro's _astro/);
  * it serves only web media types, so the gzip'd .bin.gz data become base64 text
    (*.b64.txt), which LandfallFigure decodes;
  * the unused demo folders and unused data cases are dropped.
Usage:  python3 scripts/make-artifact.py OUT_DIR      (then publish OUT_DIR/index.html)
        PALETTE=mono python3 scripts/make-artifact.py OUT_DIR   (monochrome comparison copy)
Prints the list of supporting files to pass to the Artifact publish call.
The build goes to OUT_DIR.build, so dist/ and dist-mono/ are left alone.
"""

import base64
import glob
import json
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(out):
    env = dict(os.environ, SITE_BASE="/__B__")  # PALETTE passes through to the build
    build = out + ".build"
    shutil.rmtree(build, ignore_errors=True)
    subprocess.run(["npx", "astro", "build", "--outDir", build], cwd=ROOT, env=env, check=True, capture_output=True)
    shutil.rmtree(out, ignore_errors=True)
    shutil.copytree(build, out)
    shutil.rmtree(build)
    os.chdir(out)
    for dp, _, fs in os.walk("."):  # never publish macOS metadata
        for f in fs:
            if f == ".DS_Store":
                os.remove(os.path.join(dp, f))
    for d in ("animations", "figures", "videos"):
        shutil.rmtree(d, ignore_errors=True)
    os.rename("_astro", "assets")

    h = open("index.html").read().replace("/__B__/_astro/", "./assets/").replace("/__B__/", "./")
    head = re.search(r"<head>(.*)</head>", h, re.S).group(1)
    body = re.search(r"<body[^>]*>(.*)</body>", h, re.S).group(1)
    title = re.search(r"<title>.*?</title>", head, re.S).group(0)
    head = re.sub(r'<meta charset="[^"]*"\s*/?>|<meta name="viewport"[^>]*>', "", head.replace(title, ""))
    # the host supplies <html>, so carry the dev-only palette switch over in a script
    palette = re.search(r'<html[^>]*\bdata-palette="([^"]+)"', h)
    if palette:
        head = f'<script>document.documentElement.dataset.palette = "{palette.group(1)}";</script>\n' + head
    open("index.html", "w").write(title + "\n" + head.strip() + "\n" + body.strip() + "\n")
    for css in glob.glob("assets/*.css"):
        s = open(css).read()
        open(css, "w").write(s.replace("/__B__/_astro/", "./"))
    # legacy KaTeX font formats (every current browser uses .woff2); keep fonts that lack a .woff2
    for f in glob.glob("assets/KaTeX_*.ttf") + glob.glob("assets/KaTeX_*.woff"):
        if not glob.glob(f.rsplit(".", 2)[0] + ".*.woff2"):
            continue
        os.remove(f)

    # keep only the data folders the page uses (data-base="./data/<case>/"), converted to base64 text
    used = set(re.findall(r'data-base="\./data/([^/"]+)/"', h))
    for case in os.listdir("data"):
        if case not in used and os.path.isdir(os.path.join("data", case)):
            shutil.rmtree(os.path.join("data", case))
    for case in used:
        data = os.path.join("data", case)
        man_path = os.path.join(data, "manifest.json")
        man = open(man_path).read()
        for f in glob.glob(os.path.join(data, "*.bin.gz")):
            name = os.path.basename(f)
            new_name = name.replace(".bin.gz", ".b64.txt")
            open(os.path.join(data, new_name), "w").write(base64.b64encode(open(f, "rb").read()).decode())
            os.remove(f)
            man = man.replace(f'"{name}"', f'"{new_name}"')
        open(man_path, "w").write(man)

    leftovers = subprocess.run(["grep", "-rl", "__B__", "."], capture_output=True, text=True).stdout
    if leftovers:
        sys.exit(f"unresolved base paths in: {leftovers}")
    files = sorted(p[2:] for p in (os.path.join(dp, f) for dp, _, fs in os.walk(".") for f in fs) if p != "./index.html")
    print(json.dumps([{"path": p} for p in files]))


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]))
