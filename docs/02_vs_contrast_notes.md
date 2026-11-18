# How much should Vs drop at the BDT?

Notes behind `VS_RATIO` (Vs drop / Vp drop) in `04_scenarios` and `05_matrix_imaging`, used by `05b_contrast_sweep`.

## Background velocities

The model uses the Sewell et al. (2017) 1D model (Mroczek et al. 2020, Table 2). Near the BDT depths (4–6 km) it
gives Vp ≈ 5.4 km/s, Vs ≈ 3.2 km/s and Vp/Vs ≈ 1.69. That fits the regional tomography: Sherburn et al. (2003)
found Vp/Vs falling from > 1.76 to 1.70–1.73 between 4 and 6 km over much of the central TVZ, and noted that Vp/Vs
anomalies are much smaller than Vp anomalies. Bench measurements on Torlesse greywacke (sealed sandstone, room
conditions) give Vp 5.26 and Vs 3.17 km/s, so Vp/Vs ≈ 1.66.

Nothing published measures a velocity step *at* the BDT under Rotokawa. The Vs drop has to come from what is likely
to cause the drop, so it is tied to the Vp drop through the ratio below rather than set on its own.

## Ratio for each cause

| Cause | Vp drop | Vs drop | Vs/Vp drop | Vp/Vs | Basis |
|---|---|---|---|---|---|
| Temperature alone, per 100 °C | 0.3–0.9 % | 0.4–1.2 % | 1.3–1.5 | slightly up | Kern & Richter (1981): dVp/dT −1.6 to −4.9 and dVs/dT −1.4 to −3.9 ×10⁻⁴ km/s/°C, 20–500 °C at 6 kbar, scaled by Vp 5.4 and Vs 3.2 km/s |
| Fluid-filled thin cracks (liquid-like brine or dense supercritical fluid) | 0.20 ε | 0.46 ε | **2.25** | up | O'Connell & Budiansky (1974), dilute limit, Poisson's ratio 0.25, crack density ε |
| Partial melt | — | — | ≥ 2 | up | same physics as fluid-filled inclusions; Bannister et al. (2004) read low-Vs layers at depth in the TVZ as melt |
| Vapour-like fluid (or empty) cracks | 1.25 ε | 0.72 ε | **0.58** | down | O'Connell & Budiansky (1974), dry cracks, same assumptions |

Crack derivation (dilute, ν = 0.25, so K = 5/3 G):

- dry: ΔK/K = −16/9 (1−ν²)/(1−2ν) ε = −3.33 ε, ΔG/G = −32/45 (1−ν)(5−ν)/(2−ν) ε = −1.45 ε
- saturated: ΔK/K ≈ 0, ΔG/G = −32/15 (1−ν)/(2−ν) ε = −0.91 ε
- ΔVp/Vp = ½ ΔM/M with M = K + 4G/3; ΔVs/Vs = ½ ΔG/G (density change from thin cracks neglected)

A crack density of 0.05–0.1 gives Vp drops of 1–2 % and Vs drops of 2.3–4.6 % for fluid-filled cracks.

## Which applies at Rotokawa

- **Temperature alone is small.** Even 100 °C across the top of the BDT gives only about 1 % in Vs. It sets a
  floor, not the likely contrast.
- **Fluid in cracks is the likely cause, and it should act like a liquid.** At 4–6 km the pore pressure is at
  least hydrostatic (~40–60 MPa) and may approach lithostatic (~100–150 MPa) below the BDT. At those pressures water
  at 400–450 °C is dense (roughly 0.4–0.7 g/cm³; approximate, check against IAPWS-95). It is still stiff
  compared with the aspect ratio of thin cracks, so the cracks behave as saturated: Vs drops about twice as much as
  Vp, and Vp/Vs rises.
- **The vapour-like case belongs to shallower, low-pressure zones,** such as the steam zones at Larderello, or the
  Vp/Vs < 1.4 zone at 2–3 km under Kuju. It is kept in the sweep as the end-member where P is the better phase.
- **Melt**, if the deep conductor (Heise et al. 2008) is partly molten, also gives a Vs drop of at least twice the
  Vp drop.

**Default:** `VS_RATIO = 2.0` (a little under the 2.25 crack value, to allow some temperature contribution). With
`DROP_VP = 0.05` this gives a 10 % Vs drop and raises Vp/Vs inside the zone from ~1.69 to ~1.78.

## Caveats

- These are dilute, non-interacting crack formulas; at ε > ~0.1 they overestimate the drop.
- The sign is not certain. If cracks close below the BDT and pore pressure stays low, velocities could rise instead.
  The LVZ model assumes overpressured fluid keeps cracks open.
- S imaging in 05 is the same scalar Born model with Vs in place of Vp. It ignores P–S conversions, density changes
  and S-wave polarisation, and assumes the ambient-noise correlations give S body waves as cleanly as P.

## References

- Bannister, S., Bryan, C. J., & Bibby, H. M. (2004). Shear wave velocity variation across the Taupo Volcanic Zone,
  New Zealand, from receiver function inversion. *Geophysical Journal International*, 159, 291–310.
- Heise, W., Caldwell, T. G., Bibby, H. M., & Bannister, S. C. (2008). *Geophysical Journal International*, 173(2), 740–750.
- Kern, H., & Richter, A. (1981). Temperature derivatives of compressional and shear wave velocities in crustal and
  mantle rocks at 6 kbar confining pressure. *Journal of Geophysics*, 49, 47–56.
- O'Connell, R. J., & Budiansky, B. (1974). Seismic velocities in dry and saturated cracked solids. *Journal of
  Geophysical Research*, 79(35), 5412–5426.
- Sherburn, S., Bannister, S., & Bibby, H. (2003). Seismic velocity structure of the central Taupo Volcanic Zone,
  New Zealand, from local earthquake tomography. *Journal of Volcanology and Geothermal Research*, 122, 69–88.
- Greywacke bench data: *Greywacke Fault Zone Elastic Properties: Field and Laboratory Constraints*, PhD thesis,
  Te Herenga Waka–Victoria University of Wellington.
