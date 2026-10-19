# VERIFICATION_MERRA2 — MERRA-2 KEV/VTE budgets for the 3–5 Jan 2023 landfalling AR

Produced on Sherlock, 2026-10-13, in `$GROUP_SCRATCH/Aidi/VKE_budget/case2023/`.

Deliverable: `budget_merra2_2023-01.nc` (11.5 MB, sha256 `6e6c2e6d…66f5a0`).

## 1. Summary

- **The Python port reproduces the original NCL to float32 rounding once its four deliberate departures are switched back.** I switched them on the NCL side (§3, Method). With all four on, the NCL matches the port to a relative RMS of 1e-7 for every term; IPEKE reaches 2e-5, which is float32 cancellation in ∇Φ. With all four off, the NCL variant is bit-identical to `ncl/ARanalysisMERRA2_IVKE1.ncl`. No unexplained differences remain, so **no bugs were found in the port's translation**.
- **The departures do move the terms, by 2–30% RMS.** The vertical-integration change dominates for the boundary-layer and vertical terms. The finite-difference and isobaric changes dominate for the horizontal terms (§3, Table b).
- **The integration change fixes an NCL problem.** When surface pressure exceeds 1000 hPa (most of the ocean), `vibeta` gives the 1000 hPa level zero weight, so the original NCL drops the lowest 12.5 hPa from every column integral (§4).
- **Budget closure inside `mask_IKEV`, over all 16 frames: RMS(KEV_resid)/RMS(KEV_tend) = 0.30. For VTE inside `mask_IVTE` it is 0.32.** The residual doesn't come from the port: the original NCL closes no better (0.26 vs 0.24 on 2023-01-04). With MERRA-2's own dynamics tendencies (DUDTDYN, DQVDTDYN) the same budget closes to 3%. The residual therefore comes from recomputing the dynamics offline from 3-hourly, 0.5° output (§6).
- **KEV_ana** has an RMS of 0.038 kg s⁻³ inside the mask, about 11% of RMS(KEV_tend) and of RMS(KEV_hflux). Its mean is about 0 and its local maximum is 0.31 kg s⁻³ against a tendency maximum of 2.3. My recommendation: don't give it a headline panel; show it as a minor term or fold it into an "other" group (§6).
- **VTE budget:** implemented with your exact HCV/VCV form. The NCL's commented-out form changes the horizontal + vertical sum by 6.7% RMS, not zero, because the two forms differ by q²K∇₃·v, which vanishes only by continuity (§5).
- **Ong & Yang's `ARanalysisMERRA2.ncl` isn't an independent check.** It's Zhang's NCL with the q²K lines switched on, using the same `vibeta`/`gradsf` code. It reproduces my uncommented NCL VTE variant bit for bit (§5).

## 2. What was run

| step | script | notes |
|---|---|---|
| download | `case2023/download_MERRA2_2023.sbatch` | 13 files, 16 GB, into `case2023/data/2023` (`../data` belongs to Aidi and isn't group-writable); HDF5 magic and size check: 13/13 ok |
| environment | `case2023/env_vke.sh` | conda env `$GROUP_HOME/yucli/envs/vke`: Python 3.12, numpy 2.0.0, xarray 2025.1.1, metpy 1.6.3, tempest-extremes, ffmpeg; NCL from module `devel ncl/6.6.2`. Aidi's venv isn't readable by others. |
| KEV + VTE budget | `case2023/py/ARanalysisMERRA2_IVKE1_VTE.py --vte --complevel 4` via `case2023/run_port.sbatch` | a copy of `py/ARanalysisMERRA2_IVKE1.py`, importing the original helpers from `../py`. Without `--vte` its output is bit-identical to the original (checked on 2017-07-12: all 38 2D fields and the inst file). About 5.5 min per day, 40 GB request, 4.3 GB tavg file per day. |
| AR masks | DetectBlobs in `run_port.sbatch`, same command as stage 4 of `run_example.sbatch` | `_LAPLACIAN{9,7.0}(X) < thr`, ≥ 50 points, \|lat\| ≥ 20°; thresholds from `threshold_v2.txt`: Lap(IVT) −41023.514894, Lap(IKEV) −478203.362115, Lap(IVTE) −2105.315636. My conda build reproduces Aidi's 2017-07-12 masks bit for bit. No TC filter. |
| website file | `case2023/py/make_website_budget.py` | see §7 |
| verification | `case2023/verify_ncl.sbatch` (NCL, 64 GB), `case2023/ncl/ARanalysis_verify.ncl`, `case2023/verify/{compare,report_tables}.py` | full tables in `case2023/verify/tables_20230104.md` and `tables_20170712.md` |

Paths: none of the new scripts point into Aidi's tree. They use `$GROUP_SCRATCH/Aidi/VKE_budget` and run under account `dayang` on partition `serc`. Aidi's own files were not modified.

## 3. Verification: Python port vs original NCL (2023-01-04; 2017-07-12 repeated in the appendix file)

### Method

The port can't be switched back to NCL behavior: it has no `vibeta` and no spherical-harmonic package. So I worked from the NCL side. `case2023/ncl/ARanalysis_verify.ncl` is the original NCL algorithm with three switches; the FD switch carries two of the port's departures, the finite differences and the `divergence_2d` metric term:

- `integ`: `vibeta` vs the port's trapezoid. The trapezoid starts at the lowest above-ground level and drops any layer with a missing end point, exactly as `tools_comp.integZ` does.
- `hgrad`: `gradsf`/`uv2dvF`/`grad_latlon_cfd` vs the port's finite differences. The port's FD is metpy `geospatial_gradient`. Its per-latitude weights are exported from metpy itself, and they reproduce `tools_comp.gradient_2d` and `divergence_2d` to 1e-11.
- `order`: horizontal derivatives taken on model levels and then interpolated, vs taken on pressure levels.

It also computes the VTE terms. Configuration tags: `v000` = the NCL algorithm, `v111` = the port algorithm. The original script (`orig`) and Ong & Yang's script (`ong`) were run unmodified.

Statistics are over |lat| < 87° (the port blanks poleward), all 8 tavg steps, with B as the reference. The tables below are global; the website-domain versions (190–260°E, 15–60°N) are in `verify/tables_20230104.md` and give the same picture.

**(a) Port vs original NCL**

| field | max abs diff | rel RMS | corr |
|---|---|---|---|
| IKEV (inst) | 1.72e+03 | 3.76e-02 | 0.999107 |
| IVT (inst) | 1.22e+02 | 5.01e-02 | 0.997781 |
| IVTE (inst) | 3.11e+01 | 6.43e-02 | 0.997531 |
| IKEVTEND / IVKE1TEND | 1.62e-01 | 7.69e-02 | 0.997165 |
| I_KEVTEND / I_VKE1TEND | 1.42e-01 | 4.42e-02 | 0.999145 |
| SPTE | 7.45e-09 | 6.92e-08 | 1.000000 |
| IHAKE | 1.41e+00 | 2.80e-01 | 0.964367 |
| IVAKE | 2.13e-01 | 1.29e-01 | 0.991990 |
| IPEKE | 8.34e-01 | 1.93e-01 | 0.982516 |
| IHAV | 1.49e+00 | 2.46e-01 | 0.969537 |
| IVAV | 7.93e-02 | 1.37e-02 | 0.999905 |
| IHCV | 1.44e+00 | 2.31e-01 | 0.973624 |
| IVCV | 2.39e-01 | 1.50e-01 | 0.989802 |
| IMOKE | 2.74e-02 | 2.90e-01 | 0.957364 |
| ITOKE | 4.13e-01 | 2.64e-01 | 0.965582 |
| IGOKE | 2.62e-02 | 1.13e-01 | 0.993755 |
| IPOV | 7.99e-02 | 1.10e-02 | 0.999938 |
| ITOV | 1.86e-02 | 7.88e-02 | 0.996545 |
| ICOV | 1.36e-12 | 6.57e-08 | 1.000000 |
| IAOKEV / IAOVKE | 3.65e-02 | 6.44e-02 | 0.998019 |

Units: kg s⁻² for IKEV and IVTE, kg m⁻¹ s⁻¹ for IVT, kg s⁻³ for the rest. ICOV is zero in MERRA-2 for this period.

**(b) Change made by each switch**, as relative RMS: RMS(A−B)/RMS(B)

| field | v000 vs orig | integration v100−v000 | FD + WGS84 v110−v100 | isobaric v111−v110 | port − v111 |
|---|---|---|---|---|---|
| IKEV (inst) | 0 | 3.8e-02 | 0 | 0 | 9.9e-08 |
| IKEVTEND | 0 | 7.7e-02 | 0 | 0 | 4.2e-07 |
| I_KEVTEND | 0 | 4.4e-02 | 0 | 0 | 1.0e-07 |
| SPTE | 0 | 0 | 0 | 0 | 6.9e-08 |
| IHAKE | 0 | 5.9e-02 | 2.1e-01 | 1.5e-01 | 2.8e-07 |
| IVAKE | 0 | 1.3e-01 | 0 | 0 | 1.0e-07 |
| IPEKE | 0 | 1.9e-01 | 2.3e-02 | 0 | 2.2e-05 |
| IHAV | 0 | 1.9e-02 | 1.4e-01 | 1.8e-01 | 3.0e-07 |
| IVAV | 0 | 1.4e-02 | 0 | 0 | 9.0e-08 |
| IHCV | 0 | 6.2e-02 | 1.8e-01 | 1.3e-01 | 2.9e-07 |
| IVCV | 0 | 1.5e-01 | 0 | 0 | 1.1e-07 |
| IMOKE | 0 | 2.9e-01 | 0 | 0 | 1.5e-07 |
| ITOKE | 0 | 2.6e-01 | 0 | 0 | 1.2e-07 |
| IGOKE | 0 | 1.1e-01 | 0 | 0 | 1.0e-07 |
| IPOV | 0 | 1.1e-02 | 0 | 0 | 8.8e-08 |
| ITOV | 0 | 7.9e-02 | 0 | 0 | 1.1e-07 |
| IAOKEV | 0 | 6.4e-02 | 0 | 0 | 1.1e-07 |

How to read it:

- **v000 vs orig = 0:** my NCL variant reproduces the original exactly.
- **port − v111 ≈ 1e-7:** the port is the switched NCL, to float32 rounding.
- **Each switch only touches the terms it should.** The FD and isobaric switches affect only the horizontal-derivative terms; the integration switch affects all column integrals. The one-at-a-time runs (`v010` = FD alone, `v011` = isobaric with `vibeta`) agree with the cumulative path.
- **The `divergence_2d` metric fix** isn't a separate NCL switch. The NCL's `uv2dvF` is spherical-harmonic and so already contains the V·tanφ/a term; the fix only matters relative to Zenodo's `tools_comp.py`, which omitted the term. Its effect is inside the FD + WGS84 column for IHCV.

**Why the integration switch is large for IMOKE, ITOKE, IVAKE and IVCV (13–29%):** those integrands peak at the lowest levels, which is exactly where the NCL's `vibeta` call drops the 1000–987.5 hPa half-layer when ps > 1000 hPa (§4.1). It's a correction, not just a difference. The column fields (IKEV, IVT, IVTE) move by 4–6% for the same reason.

The 2017-07-12 run reproduces all of the above (`verify/tables_20170712.md`): v000 vs orig = 0, port − v111 ≤ 3e-7 (IPEKE 3e-5), with switch effects of the same sizes.

## 4. Things in the code that look wrong, not just different

### Original NCL (`ncl/ARanalysisMERRA2_IVKE1.ncl`; also Ong & Yang's script)

1. **`vibeta` drops the bottom half-layer whenever ps > pbot.** In a one-column test with the script's 24 levels, pbot = 1000 hPa and ptop = 200 hPa:
   - ∫1 dp with ps = 1013 hPa gives 78750 Pa instead of 80000 Pa;
   - with ps = 1000 hPa and with ps = 980 hPa the result is exact;
   - an integrand that is non-zero only at 1000 hPa integrates to 0.

   So over most of the ocean every column integral omits the 1000–987.5 hPa half-layer (test in `verify/vibeta_test/t.ncl`). Biggest effect: surface-drag and turbulence terms.
2. **IPEKE's pole rows hold IVKE1AVG** (it's initialized with `IPEKE = IVKE1AVG` and only lat(1:n−2) is overwritten): a different quantity, in different units.
3. **`vibeta` aborts (ier = −999) on any all-missing column**, which is why IPEKE is integrated over lat(1:n−2) only. The pole rows of PEKE are missing because of `grad_latlon_cfd`.
4. **Mixed horizontal operators.** HAKE, HAV and HCV use spherical harmonics on terrain-following model levels; PEKE uses cyclic centred differences of Φ on isobaric surfaces. So the advection terms and the pressure-work term of the same budget are differenced on different surfaces, with different schemes.
5. The VTE lines that are commented out use HCV = −2qK∇·(qu) (§5). That form is internally consistent, but it isn't −∇·(u q²K).
6. Cosmetic: the `long_name` of `I_VKE1TEND` reads "Integral of VKE tendency", which refers to q²K, not q·K.

### Python port

7. **metpy double-counts the ellipsoid meridionally.** `tools_comp.gradient_2d` calls metpy `geospatial_gradient` with no CRS, which falls back to `+proj=latlon` (WGS84). metpy takes dy as the geodesic arc length (M·dφ) and then multiplies by `meridional_scale` = a/M, so every ∂/∂y is too large by a/M: +0.67% at the equator, +0.17% at 45°. The zonal direction is consistent (factor 1/(N cosφ)).
8. **`divergence_2d` mixes Earth models:** those WGS84 derivatives plus a spherical metric term, −V tanφ / `earth_avg_radius`. Small, and probably what SETUP.md's "agree to 0.3%" refers to.
9. **The port's column bottom is the lowest above-ground pressure level.**
   - It omits the partial layer between ps and that level when ps < 1000 hPa.
   - The bottom jumps discretely when ps crosses a level, so IKEVTEND contains small jumps.
   - SPTE, the continuous Leibniz term, is therefore only approximately the right boundary term for the port's column. Numerically minor: RMS(SPTE) is 2% of RMS(KEV_tend) inside the mask.
10. Longitude derivatives are one-sided at the ±180° seam rather than cyclic (outside the website domain). The `'\i'` SyntaxWarnings in the `tools_comp.py` docstrings are cosmetic.

## 5. VTE (q²K) budget

**Definitions** (case2023 copy, `--vte`): momentum terms weighted by q², vapor terms by 2qK.

- IHAKE_v = ⟨−q² u·∇K⟩; IVAKE_v = ⟨−q²(uω∂u/∂p + vω∂v/∂p)⟩; IPEKE_v = ⟨−q² u·∇Φ⟩
- I{MO,TO,GO}KE_v = ⟨q²(u·DUDT* + v·DVDT*)⟩; I{PO,TO,CO}V_v = ⟨2qK·DQVDT*⟩; IAOKEV_v uses the same weights
- **IHCV_v = ⟨−K∇·(q²u)⟩ and IVCV_v = ⟨−K∂(q²ω)/∂p⟩, exactly as you specified.** With these, IHAKE_v + IHCV_v = −⟨∇·(u q²K)⟩ and IVAKE_v + IVCV_v = −⟨∂(ω q²K)/∂p⟩ analytically.
- VTE_tend = IVTETEND (from the instantaneous fields); I_VTETEND is also written. SPTE_v uses VTE at the bottom level.

**About the brief's NOTE:** the NCL's commented-out form, −2qK∇·(qu) and −2qK∂(qω)/∂p (written as `IHCV_vO`/`IVCV_vO`), differs from your form by +q²K∇·u horizontally and +q²K∂ω/∂p vertically. So the horizontal + vertical sum differs by q²K∇₃·v. That is zero only if continuity holds exactly; **in this data it is 6.7% of RMS(VTE_hflux + VTE_vflux).** The split changes, and so does the sum, though less. The website file uses your form.

**Port vs NCL (Table c in `verify/tables_20230104.md`):**

- Port vs `v111` (the NCL-side VTE implementation with all port switches on): ≤ 5e-7 for every term, 1.5e-5 for IPEKE_v.
- Port vs the original-algorithm NCL VTE (`v000`): 2–34% RMS, explained by the same three switches as for KEV.
- Your "check 1" (a copy of the NCL with the VTE lines uncommented and the HCV/VCV fix) is `v000`'s `*_v` variables.

**Ong & Yang:** `ARanalysisMERRA2.ncl` (github.com/HingOng/AR_IVKE, commit ffd9a4f) is Zhang's NCL with the q²K lines active. It passes integrands to `vibeta` as `x(time|:,lat|:,lon|:,lev|:)`, in the same lev order as `p` (1000 → 150 hPa), so the ordering is correct. Its output is **bit-identical** (rel RMS 0) to my `v000` VTE terms for every shared term. HCV/VCV are compared against my `_vO` pair.

Per-term relative RMS, port vs Ong (2023-01-04, global):

| term | rel RMS |
|---|---|
| IVTEAVG | 0.068 |
| IVKETEND | 0.12 |
| IHAKE | 0.32 |
| IVAKE | 0.19 |
| IPEKE | 0.22 |
| IHAV | 0.30 |
| IVAV | 0.021 |
| IMOKE | 0.31 |
| ITOKE | 0.28 |
| IGOKE | 0.13 |
| IPOV | 0.019 |
| ITOV | 0.11 |
| IAOVKE | 0.094 |

As for KEV, these differences are fully explained by the integration and derivative switches. They are not an independent confirmation of the port.

## 6. Budget closure and KEV_ana

**Over the 16 website frames, inside mask ≥ 0.5** (KEV: 15844 point-frames; VTE: 15324), RMS in kg s⁻³:

| term | KEV | VTE |
|---|---|---|
| tend | 0.354 | 2.70e-3 |
| hflux | 0.372 | 2.65e-3 |
| vflux | 0.025 | 2.3e-4 |
| pe | 0.246 | 1.75e-3 |
| cond | 0.196 | 1.60e-3 |
| fric | 0.132 | 9.1e-4 |
| ana | 0.038 | 2.3e-4 |
| spte | 0.007 | 6.0e-5 |
| resid | 0.106 | 8.6e-4 |
| **RMS(resid)/RMS(tend)** | **0.299** | **0.320** |

Over the whole domain the ratios are 0.325 (KEV) and 0.326 (VTE).

**Where the KEV residual comes from** (same mask and frames):

| variant | RMS(resid)/RMS(tend) |
|---|---|
| your definition (IKEVTEND, flux form, + SPTE) | 0.299 |
| without SPTE | 0.302 |
| advective form (IHAV + IVAV instead of IHCV + IVCV) | 0.268 |
| tend = I_KEVTEND | 0.293 |
| MERRA-2 dynamics (IDOKE + IDOV) instead of hflux + vflux + pe, tend = IKEVTEND | 0.156 |
| **MERRA-2 dynamics, tend = I_KEVTEND** | **0.028** |

- The physics, analysis and tendency diagnostics close to 3%.
- The residual comes from two sources:
  1. the offline explicit dynamics (hflux + vflux + pe) differ from MERRA-2's own DYN tendency by 23% RMS (finite differences of 3-h-mean states on the 0.5° pressure-level grid);
  2. IKEVTEND (instantaneous Np columns) differs from I_KEVTEND (tavg Nv-interpolated weights) by 15%.
- The flux-form/continuity term qK∇₃·v is 5% of the tendency.
- **The port doesn't make closure worse.** On 2023-01-04 the brief-definition ratio is 0.262 for the original NCL, 0.249 for trapz + FD, and 0.241 for the port; all variants give 0.027–0.028 with MERRA-2's dynamics.
- For VTE on 2023-01-04: 0.268 with your definition, 0.031 with MERRA-2's dynamics.

**KEV_ana:**
- RMS 0.038 vs tend 0.354, hflux 0.372 and pe 0.246; mean −0.0003 (no systematic sign).
- Median |ana| / Σ|terms| is 0.038; maximum |ana| is 0.31 kg s⁻³ against a maximum |tend| of 2.3.
- VTE_ana / VTE_tend is about 8%.
- **Recommendation:** keep it in the data and the residual definition, but don't give it its own headline panel. Grouping it with `spte` as "analysis + surface pressure" is honest at this magnitude.

## 7. Output file `budget_merra2_2023-01.nc`

- **Dimensions:** `time` (16) × `latitude` (91, 15–60°N) × `longitude` (113, 190–260°E, converted from MERRA-2's −180…179.375).
- **Grid and storage:** native 0.5° × 0.625°, float32, zlib level 4 with shuffle, NaN fill.
- **Time:** frame centres 2023-01-03 13:30Z … 2023-01-05 10:30Z, every 3 h, in `hours since 2023-01-01`, which decodes to datetime64. Each time stamp is the centre of a MERRA-2 tavg window. `time_bnds_inst_lo`/`_hi` give the two bounding instantaneous times (tc ± 1h30), which are also the two steps IKEVTEND and IVTETEND difference.
- **Column fields** `IVT` (kg m⁻¹ s⁻¹), `IVTE` and `IKEV` (kg s⁻²): mean of the two bounding instantaneous steps.
- **Budget fields** `KEV_{tend,hflux,vflux,pe,cond,fric,ana,spte,resid}` and `VTE_{…}` (kg s⁻³): 3-h means of 1000–200 hPa column integrals, defined as in the brief.
  - Matching the brief, `KEV_hflux = IHAKE + IHCV`.
  - The port's own flux-form variables for the same quantity are `IHFOKE + IHFOV + IHFADT` (and `IVFOKE + IVFOV + IVFADT` for vflux). Equal analytically; discretely 5.8% (h) and 7.0% (v) RMS apart, because the products are differenced in a different order. No single port variable equals it.
- **Masks** `mask_IVT`, `mask_IVTE`, `mask_IKEV`: mean of the two bounding 0/1 DetectBlobs masks, so 0, 0.5 or 1. The 0.5 contour encloses points flagged at either bounding time.
- **Attributes:** script per field (`source`), time meaning, thresholds and percentile levels, detection command, grid and integration notes, and `code_sha256_16` (hashes of the port, its helpers, the VTE copy, the assembly script and `threshold_v2.txt`; there's no VCS).
- **Size:** 11.5 MB rather than the 50–100 MB estimate, since the domain is 91 × 113 points.

## 8. Reproduce

```
cd $GROUP_SCRATCH/Aidi/VKE_budget/case2023
sbatch download_MERRA2_2023.sbatch                                   # needs ~/.netrc for urs.earthdata.nasa.gov
DATE=20230104 sbatch --export=ALL verify_ncl.sbatch                  # 8 NCL configs, 64 GB, ~20 min each
for d in 20230103 20230104 20230105; do DATE=$d OUT=output EXTRA="--complevel 4" sbatch --export=ALL run_port.sbatch; done
source env_vke.sh && python py/make_website_budget.py --stats verify/website_stats.json
cd verify && python report_tables.py 20230104 ../output > tables_20230104.md
```

The large intermediates stay on Sherlock: `case2023/output/*tavg*` (3 × 4.3 GB), `case2023/verify/ncl_*/` (7.7 GB per day) and `case2023/data/2023` (16 GB). They're all on purged scratch (90 days); copy them to `$GROUP_HOME` or Oak if they need to persist.
