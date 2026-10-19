# Verification against the authors' code

All comparisons ran on the same ERA5 input (ARCO, 0.25°, 25 pressure levels from 1000 to 150 hPa, lon 160–285°E, lat 5–72°N). Every remaining difference therefore comes from the code. Scripts are in `verify/`, and the authors' original code is in `../reference/`.

## Summary

| Check | Result |
|---|---|
| **A.** Ong & Yang `ARanalysisERA5.ncl` (run in NCL 6.6.2) vs our port | **Bit-for-bit identical:** 9 fields × 6 times, including the landfall. One field differs by 1 ulp at a few points |
| **A.** Our production budget vs the verified port | Agree to float32 rounding (relative RMS 6–9×10⁻⁸) for every shared term |
| **B.** Zhang et al. KEV lines (NCL, global ERA5) vs ours | Identical formulas agree to rounding. Spectral vs finite-difference derivatives differ by 8–10% (corr. 0.996) |
| **C.** TempestExtremes `DetectBlobs` vs our port | **Identical masks:** 0 of 59,051 AR cells differ |
| **D.** Ong & Yang Fig. 5 growth rates (MERRA-2) vs ours (ERA5) | 3 Jan: +0.76/day vs +0.91. 4 Jan: −0.43/day vs −0.75 (see D) |
| NCL primitives (`vibeta`, `center_finite_diff_n`, `grad_latlon_cfd`) | Bit-exact against NCL on 4,000 random columns and a full test grid |

## Discrepancies found in the authors' code

1. **Vertical integrands are reversed in `ARanalysisERA5.ncl`. This affects Ong & Yang's ERA5 results only.**
   - **The mechanism:**
     - The script reads `level({1000:150})`. On NCAR RDA ERA5 files, `level` is stored ascending (1 → 1000 hPa; checked on the actual file `e5.oper.an.pl.128_133_q.ll025sc.2023010300_2023010323.nc`). NCL therefore returns 1000 … 150, so `p` runs bottom-to-top, which is the order `vibeta` requires.
     - Every integral except `IVAV` passes its integrand as `X(...,level|::-1)`, which runs top-to-bottom. The integrand is therefore paired with the wrong pressures: the value at 150 hPa gets the 1000 hPa layer weight, and so on.
     - Over terrain, the surface mask also removes the wrong levels. Extrapolated below-ground values get integrated, and the upper levels are dropped instead.
   - **Effect on this case** (faithful vs. corrected): IVT 23%, IVKE (= IVTE) 29–30%, IHAKE 58–60%, IVAKE 61–69%, IPEKE 25–28%, IHAV 47–52%, IKEpart 46–50% and IVPpart 47–49% relative RMS. Spatial correlation stays at 0.90–0.99.
   - **Scope:** this affects the ERA5 robustness checks (Fig. 4c,d, Supplementary Fig. S9), and any ERA5-based detection run on these fields. Both MERRA-2 scripts (Ong & Yang's and Zhang et al.'s) pass integrands in the correct order (`lev|:`), so the papers' main results are **not** affected.
   - **What we do:** `ong_compat.analyze(faithful=True)` reproduces the published behavior, and the production code uses the correct order.

2. **`ARanalysisERA5.ncl` cannot run as published.** It reads `Q(0:n_tim-1)` but indexes `U(2:n_tim+1)` and `Q(1:n_tim)` for the centred tendencies, which is out of range. The smallest fix is four marked edits, saved as `reference/ong_yang_2024/ARanalysisERA5.patch`: carry one extra hour on each side, and read those extra hours separately for the tendency terms. Nothing else was changed.

3. **Time alignment in `ARanalysisERA5.ncl`.** With that fix, the tendency parts (`IKEpart`, `IVPpart`) are centred one hour after the instantaneous fields in the same output slot. They are also integrated with the surface pressure of the earlier time. The faithful port reproduces this. The production code centres everything on the same time.

4. **Zhang et al.'s VTE horizontal/vertical split** (in the commented-out lines of `ARanalysisMERRA2_IVKE1.ncl`). It weights the vapor-flux divergence by 2qK (`HCV = -KQ2*divQ`). The result, HAKE + HCV = −q²u·∇K − 2qK∇·(qu), differs from −∇·(u q²K) by −q²K∇·u. The vertical term differs by the opposite amount, so the total transport is unchanged. For KEV their split equals the exact flux divergence. These lines are commented out, so which split was used for the published VTE figures can't be determined.

5. **Minor issue in Zhang et al.'s `AR_percentile_v3.py`.** It computes `IVKE1_std = np.std(IVKE2)` and `IVKE2_mean = np.mean(IVKE1)`, with the names swapped. These two values are printed but not used for the thresholds.

## A. Ong & Yang ERA5 script

`verify/ong_era5.py` writes an input file in RDA layout and runs the patched NCL script. It then runs `arbudget.ong_compat` in faithful and in corrected mode.

| Window | Fields | NCL vs faithful port |
|---|---|---|
| 3 Jan 11–14Z (4 output times) | IVT, IVKE, IHAKE, IVAKE, IPEKE, IHAV, IVAV, IKEpart, IVPpart | max \|Δ\| = 0 for all |
| 4 Jan 17–18Z (2 output times, landfall) | same | max \|Δ\| = 0 except IHAV (2.8×10⁻¹⁴, 1 ulp at a few points) |

To get bit-exact results, the port reproduces several details of NCL:
- `vibeta`, including its quirk that the 1000 hPa level carries zero weight whenever surface pressure is above 1000 hPa;
- the edge formulas of cyclic `center_finite_diff_n`: (q₁−q_{n−1})/2(r₁−r₀) at the start and (q₀−q_{n−2})/2(r_{n−1}−r_{n−2}) at the end;
- `grad_latlon_cfd` building its coordinates in the precision of `lat` (float32 here), including a cosine computed in double precision and then rounded.

**Production vs the verified port** (`verify/production_vs_compat.py`, VTE at 12–14Z):

| Production term | Verified-port expression | Relative RMS |
|---|---|---|
| IVTE, IVT | IVKE, IVT | 6.3×10⁻⁸ |
| pe | IPEKE | 6.4×10⁻⁸ |
| fric | IKEpart − (IHAKE + IVAKE + IPEKE) | 7.8×10⁻⁸ |
| cond | IVPpart − (IHAV + IVAV) | 9.2×10⁻⁸ |
| hflux + vflux | IHAKE + IVAKE + IHAV + IVAV | 7.9×10⁻⁸ |

These links tie each production term to fields that NCL reproduces exactly. They required two design choices:
- **Ong's discrete forms.** ∂K/∂t is computed as **u**·∂**u**/∂t, and ω∂K/∂p as ω**u**·∂**u**/∂p. Differencing K directly gives friction values about 10% different.
- **vflux defined through continuity.** It is the advective total minus the exact horizontal flux divergence. Computing −⟨∂(ωX)/∂p⟩ directly leaves a continuity defect, −⟨X(∇·u + ∂ω/∂p)⟩, equal to 31–37% of the transport. ERA5's pressure-level ω is not continuity-consistent with a 0.25° finite-difference divergence, and neither residual term would absorb that defect.

## B. Zhang et al. KEV terms

`verify/zhang_kev_era5.ncl` copies the computational lines of `ARanalysisMERRA2_IVKE1.ncl` verbatim. It changes only the inputs:
- ERA5 is already on pressure levels, so the model-level interpolation is skipped;
- ERA5's `Z` is already Φ, so it is used directly instead of `PHI = H*g`;
- latitude is ordered S → N, as `gradsf`/`uv2dvF` require.

It runs on one global hour (2023-01-04 00Z), and the comparison covers the display domain.

| Ours | Theirs | Relative RMS | Corr. | Why |
|---|---|---|---|---|
| IKEV | IVKE1 | 6.3×10⁻⁸ | 1.000 | Same formula and integral |
| pe | IPEKE | 4.7×10⁻⁵ | 1.000 | Same `grad_latlon_cfd`; latitude order changes float rounding |
| hflux | IHAKE + IHCV | 9.9×10⁻² | 0.996 | Our 0.25° regional finite differences vs their global spherical harmonics (`gradsf`, `uv2dvF`) |
| vflux | IVAKE + IVCV | 1.8 | 0.47 | Different by design (see A): theirs is the direct, near-zero column flux; ours absorbs the continuity defect |
| hflux + vflux | sum of the four | 8.1×10⁻² | 0.997 | Spectral vs finite differences, plus the continuity defect |

Terms that need MERRA-2 physics tendencies (MOKE, TOKE, GOKE, POV, TOV, COV and the analysis term) cannot be computed from ERA5. As in Ong & Yang's ERA5 analysis, we estimate them as the residuals `cond` and `fric`.

## C. Detection

`verify/detect_te.py` runs TempestExtremes `DetectBlobs` from conda-forge with Ong & Yang's exact commands (`_LAPLACIAN{9,7.0}`, thresholds −4×10⁴ for IVT and −2×10³ for IVKE, `--minabslat 20`, a 50-point filter). It adds `--regional`, because our subset is not global.

| Field | TE AR cells | Port AR cells | Mismatched |
|---|---|---|---|
| IVT (4 times) | 30,362 | 30,362 | 0 |
| IVKE (4 times) | 28,689 | 28,689 | 0 |

The port follows `BuildLaplacianOperator`:
- ring points start eastward and rotate toward north;
- each ring point snaps to its nearest grid node, with duplicates and the centre node dropped;
- each neighbour gets weight 4/(N·d²), where d is the true great-circle distance to the snapped node;
- weights are stored in float32;
- blobs use 4-neighbour connectivity.

The **IKEV threshold** is −5.5×10⁵ kg s⁻² rad⁻², the fixed reference value in Zhang et al.'s `AR_percentile_v3.py`. Their published detection used the matching climatological percentile, which needs a multi-year climatology. Over the 49 frames, the Jaccard overlap is 0.80 between IVT and IVTE masks, and 0.76 between IVT and IKEV masks.

Ong's 50-point filter was calibrated on MERRA-2's 0.5°×0.625° grid. An equivalent area filter (≥ 50 MERRA-2 cells ≈ 1.5×10⁵ km²) changes 6–7% of the AR cells on the 0.25° grid: 11,950 of 184,286 for IVT, and 12,425 of 183,988 for IVTE. The figure uses the published 50-point filter.

## D. Physics against Ong & Yang Fig. 5

The growth rate is the Eq. 11a projection over 160°W–100°W, 20°–50°N, using 3-hour means over 18–21Z (trapezoidal weights on the hourly terms) and the full IVTE field. Values are per day.

| | Paper (MERRA-2) | Ours (ERA5), from ΔI over 3 h | PE→KE | cond | fric | hflux | vflux |
|---|---|---|---|---|---|---|---|
| 3 Jan 18–21Z (intensifying) | +0.91 | **+0.76** | +3.66 | −2.20 | −1.39 | +0.40 | +0.11 |
| 4 Jan 18–21Z (after landfall) | −0.75 | **−0.43** | +3.20 | −2.16 | −2.17 | +0.60 | +0.13 |

- **3 Jan:** our value is within 17% of the paper. As in their Fig. 5c,d, PE→KE is the dominant growth term, offset by condensation and friction.
- **4 Jan:** same sign, but only 58% of the paper's magnitude. The tendency here comes straight from the change in IVTE over 3 hours, independent of the budget code, which A–C verify against the authors' code. The most likely explanation is the input data: ERA5 vs MERRA-2, plus MERRA-2's 3-hourly time-averaged fields vs ERA5's hourly instantaneous fields. The difference is largest after landfall, where the AR's decay depends on each reanalysis's precipitation, surface drag and terrain. Confirming this would mean running the pipeline on MERRA-2.

## Budget closure

Inside the AR mask over all 49 frames, the residual, tend − Σ terms, is 28% (VTE) and 20% (KEV) of the tendency's RMS. The correlation between the tendency and the sum of terms is 0.96 and 0.98. This is above the 10% we aimed for. In the Fig. 5 projections the residual is small: 0.14/day against a tendency of 0.72/day.

Splitting the residual at six frames: the surface-pressure (Leibniz) term is 2–3%. Almost all the rest, 14–37%, is the **discrete product rule in time**. Over a 2-hour centred step, Δ(q²K) ≠ 2qK·Δq + q²**u**·Δ**u** when both are evaluated at the centre time, which is Ong & Yang's discretization. Removing NCL's 1000 hPa `vibeta` quirk changes the residual by less than 1%.

An exactly closing alternative exists: time-averaged weights such as (q₊² + q₋²)/2, for which a₊b₊ − a₋b₋ = ā·Δb + b̄·Δa holds exactly. Adopting it would mean departing from the authors' discretization. **This is left as a decision for the site authors.**
