# Notebooks

| Notebook | What it does | Main outputs |
|---|---|---|
| `01_velocity_model` | Builds the 3D Vp and Vs model from the known geology: layers, faults (Calibugan et al. 2022 structure) and greywacke basement. No BDT. | `velocity_model.npz` |
| `02_designs` | Station layouts and deployment schedule: Trillium sites (balanced between hypocentre location and matrix imaging), the existing Mercury stations at Rotokawa and Ngā Tamariki (from `data/RTK_NTM_Station_Coordinates.dat`), static node designs, and the phased design (wide reconnaissance grid, then dense rolling patches). Also the station bank phase 2 can be placed on. | `stations.csv`, `designs.json` |
| `03_travel_times` | P travel times from every station to every model cell (eikonal), plus 1D tables. Slow the first time. | `tt/` |
| `04_scenarios` | Possible "true" Earths: conductor positions sampled around the drill site, combined with several BDT shapes. | `scenarios.npz`, `scenarios.json` |
| `mi_common.py` | Matrix imaging shared by 05, 05b, 06 and 07: the LVZ of a scenario as Born scatterers, reflection matrices between stations that record together, confocal imaging (3D or 1D focusing), scoring, and the adaptive phase 2. One forward model and one score, so 05-07 can be compared directly in 09. | — |
| `05_matrix_imaging` | For each scenario: phase 1 picks where the BDT high is, phase 2 images it; static designs for comparison. Scores location error and coverage of the BDT top. Noise-free. | `results_05_mi.csv` |
| `05b_contrast_sweep` | Reruns 05 (noise-free) for one or two scenarios at velocity drops of 1, 2.5, 5 and 10 %, imaging with P and with S, and finds the smallest drop the best design detects. Converts the S limit to a Vp drop for each likely cause of the LVZ (Vs drop / Vp drop, `docs/02_vs_contrast_notes.md`). | `results_05b_contrast.csv` |
| `05c_phased_schedule` | Fewer node moves: the 02 phase 1 grid, then 2-4 phase 2 deployments laid out as rolling patches, interleaved sub-grids, or a coarse "zoom" deployment then rolling patches, with deployment length varied. +5 % Vp dome, with 06-style noise (so recording time matters), picks made only from data recorded by then. Test C asks whether the array images the BDT everywhere within 4 km of the drill site (the MT target area), comparing "tiled" designs (3 deployments spread over the area, no phase 1) with the phased ones over all 15 conductor positions. Reports coverage against field trips, campaign days and node-days, and shows the images (sections, picked-top and error maps, a blank-column diagnostic and where each design is reliable). Coarser `FAST` grids by default. | `results_05c_schedule.csv`, `mi_images_05c/` |


