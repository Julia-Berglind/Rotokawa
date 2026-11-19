# Meeting notes: controls on velocity (Carolyn Boulton, 18 Nov 2026)

## Summary

- **Cracks are the biggest control on velocity.** Porosity, mineralogy, fluid type and differential pressure all matter, but less.
- How cracks are oriented, which depends on the stress field, controls the fast and slow directions. That matters for near-vertical matrix-imaging rays.
- If the rock is foliated (schist), expect **7–17 % anisotropy**, slow across the foliation.
- Put more **hydrothermal alteration** into the velocity model (see Simpson 2026). Quartz cementation is likely near the BDT.
- A flat or idealised BDT is unrealistic. **Base the BDT geometry on map-view surface heat flow**: shallower where heat flow is higher.

## 1. Primary controls on velocity

| Control | Effect | Notes |
|---|---|---|
| Cracks / fractures | Largest | Density, aspect ratio and **orientation** all matter |
| Porosity | Moderate | Mostly through pore shape: thin cracks matter much more than round pores |
| Mineralogy | Moderate | Alteration can push velocity either way (see §4) |
| Fluid type | Moderate | Liquid-like vs vapour-like changes Vs/Vp (see `02_vs_contrast_notes.md`) |
| Differential pressure / stress | Moderate | Higher effective pressure closes cracks; overpressure keeps them open |

## 2. Crack orientation and stress

**Rule:** cracks open against the least principal stress, σ₃. The crack plane contains σ₁ and σ₂, and the crack normal is parallel to σ₃.

| Vertical stress | Regime | Crack orientation | Symmetry |
|---|---|---|---|
| σ₁ | Normal faulting (rift) | Vertical, striking along σ_Hmax | HTI |
| σ₂ | Strike-slip | Vertical | HTI |
| σ₃ | Reverse, or near-lithostatic pore pressure | Horizontal | VTI |


**Rotokawa:**
- The TVZ is a rift, so above the BDT σ₁ is probably vertical and σ₃ ≈ NW–SE, at right angles to the fault strike (~030°).
- So cracks are most likely vertical and strike NE, parallel to the faults.
- Waves travelling or polarised parallel to the cracks are fast; across them, slow. S-wave splitting should have its fast polarisation NE.

**Why it matters for matrix imaging (near-vertical rays):**
- Vertical cracks: vertical P is hardly affected; vertical S splits.
- Horizontal cracks: vertical P and S are both slowed as much as possible.
- **Hypothesis:** above the BDT, vertical cracks. Below it, flat-lying cracks or fabric, from near-lithostatic fluid pressure, low differential stress, or low-angle shear fabric. That would give a velocity drop for vertical rays, and so a reflector, **from crack orientation alone**, with no change in crack density.

## 3. Foliation and anisotropy

- Analogue: Torlesse greywacke on the South Island, carried through the BDT and now **Haast Schist** along the Alpine Fault.
- Foliated rock: **7–17 % anisotropy**, **slow direction across the foliation** (aligned micas).
- Horizontal foliation would make vertical rays slow, which would be a strong contrast for matrix imaging. Steep foliation would make vertical rays fast.
- **Caveat:** Rotokawa basement (Kaweka/Torlesse) is mostly low-grade greywacke, not schist. Treat foliation as an end-member scenario, not the default.

## 4. Alteration

- See **Simpson (2026)**
- **Top ~3 km:** a lot of alteration: feldspar alteration, epidote, chlorite, etc.
  - Epidote is stiff and dense, so it raises velocity. Clays and zeolites lower it.
  - **Watch for double counting:** the Sewell et al. (2017) 1D model was fitted to local data, so it already contains some alteration. Better to vary it laterally (stronger near upflow) than to shift the whole profile.
- **At BDT depths:** Carolyn thinks the greywacke is probably being **altered to quartz and heavily cemented**.
  - Quartz has a low Poisson's ratio (Vp/Vs ≈ 1.5). Silicification raises Vp, lowers Vp/Vs, and seals cracks.
  - **This works against my current LVZ assumption** (`DROP_VP = 0.05`, `VS_RATIO = 2.0`, Vp/Vs rising). It's a physical case for a velocity *increase* at the BDT.

## 5. BDT geometry

- A uniform layer or a perfect geometric shape (my current dome, block and lens) is unrealistic.
- **Use map-view surface heat flow:** the BDT should be **shallower where surface heat flow is higher**.
- How I'd implement it:
  1. Get a surface heat flux map / grid for Rotokawa.
  2. Smooth it over a few km, because surface heat is shifted sideways by lateral flow and outflows.
  3. Map it onto the 3.5–5.5 km range in **relative** terms (highest heat flow → shallowest). Don't use a conductive formula (z = kΔT/q): the surface heat is mostly carried by moving fluid, and a conductive calculation puts the BDT far too shallow.
  4. Run variants (smoothing length, lateral offset) and keep the conductor-centred shapes for comparison.

## Actions

- [ ] Write up crack-orientation and anisotropy scenarios: horizontal cracks / VTI, foliated VTI at 7–17 %, vertical cracks / HTI
- [ ] Add a **quartz-cemented (velocity increase)** case to `docs/02_vs_contrast_notes.md` and `05b_contrast_sweep`
- [ ] Add lateral alteration to `01_velocity_model`, plus a silicified zone near the BDT, without double counting the 1D model
- [ ] Find a surface heat flow map for Rotokawa and add a heat-flow-based BDT shape to `04_scenarios`
- [ ] Find stress data for Rotokawa (borehole breakouts, focal mechanisms, e.g. McNamara et al.) to confirm σ₁ is vertical and get the σ₃ direction

## Open questions

- Is there any evidence of foliation in the Rotokawa greywacke from cores or image logs?
- Below the BDT, does fluid pressure get high enough for flat-lying cracks?
- Is there shear-wave splitting data at Rotokawa that already gives a fast direction?
