# Bar chart data for the hero figure

The right half of the hero figure (`src/components/LandfallFigure.astro`) shows, for each frame, how much each VKE budget term contributes to the AR's **growth/decay** and to its **movement**. It reads one file next to the map data:

```
public/data/<case>/budget_bars.json      (currently <case> = landfall-2023-merra2)
```

The file in the repo now holds placeholder numbers (`"filler": true`). Replace it with the real one.

## Format

```jsonc
{
  "filler": false,                        // true shows a "Placeholder numbers" tag on the chart
  "times": ["2023-01-03T13:30Z", ...],    // must equal "times" in manifest.json (same frames, same order)
  "metrics": {
    "IVTE": {
      "growth": {
        "units": "% h⁻¹",                 // free text, shown in the chart heading and tooltip
        "terms": {
          "tend":  [ ...nt numbers ],     // total tendency
          "hflux": [ ... ],               // horizontal advection
          "vflux": [ ... ],               // vertical advection
          "pe":    [ ... ],               // PE → KE conversion
          "cond":  [ ... ],               // condensation only (not evaporation or chemistry)
          "fric":  [ ... ]                // turbulence only (not gravity-wave drag or convection)
        }
      },
      "movement": { "units": "...", "terms": { same six keys } }
    },
    "IKEV": { same structure }
  }
}
```

- `nt` = number of frames in `manifest.json` (16 for the MERRA-2 case).
- Each number is the term's contribution at that frame. Positive means the term **grows** the AR (growth chart) or **drives its motion** (movement chart); negative means decay or opposing the motion.
- Use `null` for frames where a number does not exist, e.g. before the AR is detected; that bar shows "n/a".
- `tend` is the total tendency, drawn above the five process terms so readers can compare it with their sum.
- Units may differ between the two charts. Each chart keeps one axis range for all frames (the largest `|value|` in the file), so bar lengths compare across time.
- No error bars for now.

## Map files

The map chips on the left read "Condensation" and "Turbulence". For consistency with the bars, regenerate `ikev_cond` / `ivte_cond` from the condensation tendency only, and `ikev_fric` / `ivte_fric` from the turbulence tendency only (same keys and encoding as before in `manifest.json`).
