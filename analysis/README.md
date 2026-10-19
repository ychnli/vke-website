# VKE budget pipeline (ERA5) for the hero figure

This pipeline produces `public/data/landfall-2023/`, the data behind the hero figure: the atmospheric river that made landfall on the US West Coast on 3–5 January 2023 (Ong & Yang 2024, Fig. 5), recomputed with hourly ERA5.

## Setup

```sh
cd analysis
uv sync                      # native arm64 Python env in .venv
```

The verification scripts also need NCL and TempestExtremes. NCL has no Apple Silicon build, so they live in a separate Intel (osx-64) conda env that runs under Rosetta:

```sh
CONDA_SUBDIR=osx-64 conda create -p ./.ref-env -c conda-forge --override-channels \
    ncl tempest-extremes xorg-libx11 xorg-libxrender xorg-libxext xorg-libxau cairo
```

`data/ncl.sh` runs NCL from that env. It closes stdin, because NCL waits on it when there is no terminal.

## Run

```sh
uv run python scripts/01_fetch.py    # stream ARCO ERA5 (~25 GB transfer, ~2 GB kept), ~30-70 min
uv run python scripts/02_budget.py   # budget terms + AR masks -> data/budget.nc
uv run python scripts/04_export.py   # web data -> ../public/data/landfall-2023/
uv run python scripts/05_verify.py   # Fig. 5 growth rates, closure, mask checks
uv run pytest                        # unit tests
```

Verification against the authors' code is described in [VERIFICATION.md](VERIFICATION.md):

```sh
uv run python verify/ong_era5.py 2023010311 4                    # A: Ong & Yang NCL script
uv run python -m verify.production_vs_compat 2023010311 4         # A: production == verified port
uv run python -m verify.detect_te data/verify_ong/<ivt file>.nc  # C: TempestExtremes
uv run python verify/fetch_global.py 2023010400                  # B: one global hour
uv run python -m verify.zhang_kev 2023010400                      # B: Zhang et al. KEV terms
```

## Layout

| Path | Contents |
|---|---|
| `arbudget/config.py` | Case window, domains, levels, and the naming map (VTE = Ong & Yang's "IVKE") |
| `arbudget/fetch.py` | ARCO streaming, with a per-(variable, hour) cache |
| `arbudget/grid.py` | Bit-exact ports of NCL `vibeta`, `center_finite_diff_n`, `grad_latlon_cfd` |
| `arbudget/budget.py` | Production budget terms for VTE and KEV |
| `arbudget/ong_compat.py` | Line-by-line port of `ARanalysisERA5.ncl`, used only for verification |
| `arbudget/detect.py` | Port of TempestExtremes `DetectBlobs` with `_LAPLACIAN{9,7.0}` |
| `arbudget/growth.py` | Projection diagnostics (Ong & Yang Eq. 11) |
| `arbudget/export.py`, `geo.py` | Web data export and Natural Earth lines |
| `verify/` | Harnesses that run the authors' code side by side with ours |
| `../reference/` | Authors' original code (MIT), plus a minimal patch to make the ERA5 script run |
