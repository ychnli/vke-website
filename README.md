# Atmospheric Rivers & Vapor Kinetic Energy: website prototype

A static website that explains the vapor kinetic energy (VKE) diagnostics for atmospheric rivers in plain language. It is based on:

- Ong & Yang (2024), *Nature Communications* 15, 9428. [doi:10.1038/s41467-024-53369-0](https://doi.org/10.1038/s41467-024-53369-0)
- Zhang, Yang, Ong & Tan (2026), *GRL* 53, e2025GL119841. [doi:10.1029/2025GL119841](https://doi.org/10.1029/2025GL119841)

It is built with [Astro](https://astro.build) and MDX. LaTeX is rendered by [KaTeX](https://katex.org) at build time, so the published pages are plain HTML and CSS and need no math JavaScript.

## Quick start

```sh
npm install        # once
npm run dev        # live-reloading preview at http://localhost:4321
npm run build      # static site → dist/
npm run preview    # serve dist/ locally
```

## Layout

```
src/pages/index.mdx            Main page content (Markdown + LaTeX + components)
src/layouts/Base.astro         Page shell (header, footer, KaTeX CSS)
src/components/
  Figure.astro                 Static image with caption
  Video.astro                  MP4/WebM animation with native controls
  FrameAnimation.astro         Image-sequence player (play/pause, slider, speed)
  Callout.astro                Highlighted "key idea" box
src/styles/global.css          Site styles (light and dark themes)
public/figures/                Static figures         → referenced as "figures/…"
public/videos/                 Video files            → "videos/…"
public/animations/<name>/      Frame sequences        → "animations/<name>/…"
scripts/make-demo-frames.mjs   Generates placeholder frames (delete when real ones exist)
```

To add a page, create `src/pages/<name>.mdx`. It will be served at `/<name>/`.

## Writing content

**Math.** Use `$...$` for inline math and `$$...$$` for display math, the same as in LaTeX:

```mdx
The budget is

$$
\frac{\partial \mathcal{K}}{\partial t} = \underbrace{\mathcal{A}}_{\text{advection}} + \cdots
$$
```

KaTeX supports most of the math commands in amsmath. See the [list of supported functions](https://katex.org/docs/supported). If an expression fails to render, the build shows it in red rather than failing.

**Figures.**

```mdx
<Figure src="figures/vke_climatology.png" alt="Describe the figure for screen readers" label="Figure 1">
  Caption, which can include $\mathrm{J\,m^{-2}}$.
</Figure>
```

**Video animations.** Supply WebM and/or MP4 files; the browser plays the first format it supports.

```mdx
<Video sources={['videos/ar_event.webm', 'videos/ar_event.mp4']} poster="videos/ar_event.jpg" label="Movie 1">
  Caption.
</Video>
```

Add `loop autoplay` to get silent, GIF-like playback. Autoplay is turned off automatically for readers who have reduced motion enabled.

The settings below produce a web-friendly MP4 from matplotlib frames:

```sh
ffmpeg -framerate 8 -i frame_%03d.png -c:v libx264 -pix_fmt yuv420p -crf 23 -movflags +faststart ar_event.mp4
```

**Frame-sequence animations.** Put numbered images in `public/animations/<name>/`:

```mdx
<FrameAnimation
  pattern="animations/ar_event/frame_{i}.png" count={48} pad={3} fps={6}
  labels={['2017-02-06 00Z', /* one label per frame */]}
  alt="What the animation shows" label="Animation 1">
  Caption.
</FrameAnimation>
```

You can pass `frames={[...]}` with explicit paths instead of `pattern`. When the player has keyboard focus, Space plays or pauses and the ← and → keys step one frame. To keep pages fast, save frames as PNG or WebP at about 1200 px wide.

## Deploying

`npm run build` writes a self-contained static site to `dist/`. Upload the contents of `dist/` to any web server.

If the site will be served from a sub-directory (for example `https://web.stanford.edu/~user/ar/`), set the base path when you build:

```sh
SITE_URL=https://web.stanford.edu SITE_BASE=/~user/ar npm run build
```

All components resolve asset paths with `src/lib/asset.ts`, so links keep working under any base path. If you write a raw `<img>` or `<a>` that points into `public/`, use `asset()` for it as well.

### GitHub Pages

`.github/workflows/pages.yml` builds and publishes the site on every push to `main`. It reads the URL and base path from the repository, so nothing needs to be set by hand. It publishes the monochrome palette; change `PALETTE: mono` to `PALETTE: color` near the top of the workflow to publish the color one. One-time setup: push the repository to GitHub, then in **Settings → Pages** set **Source** to **GitHub Actions**. Runs show up under the **Actions** tab; **Run workflow** there redeploys without a new commit.

To check a Pages build locally (replace `<user>` and `<repo>`):

```sh
SITE_URL=https://<user>.github.io SITE_BASE=/<repo> npm run build:mono
npm run preview:mono   # then open http://localhost:4322/<repo>/
```

### Claude artifact

The artifact host needs a differently packaged copy (relative paths, no `_astro/` folder, data as base64 text); `scripts/make-artifact.py` makes it.

```sh
npm run artifact        # color site → artifact/
npm run artifact:mono   # monochrome comparison copy → artifact-mono/
```

The script prints the list of supporting files. Publishing is done from Claude: ask it to publish `artifact/index.html` (or `artifact-mono/index.html`) with those files to the existing artifact link.
