# Notebooks

| Notebook | What it does | Main outputs |
|---|---|---|
| `01_velocity_model` | Builds the 3D Vp and Vs model from the known geology: layers, faults (Calibugan et al. 2022 structure) and greywacke basement. No BDT. | `velocity_model.npz` |
| `02_designs` | Station layouts and deployment schedule: Trillium rings, the existing Mercury stations at Rotokawa and Ngā Tamariki (from `data/RTK_NTM_Station_Coordinates.dat`), static node designs, and the phased design (wide reconnaissance grid, then dense rolling patches). Also the station bank phase 2 can be placed on. | `stations.csv`, `designs.json` |
| `03_travel_times` | P travel times from every station to every model cell (eikonal), plus 1D tables. Slow the first time. | `tt/` |
| `04_scenarios` | Possible "true" Earths: conductor positions sampled around the drill site, combined with several BDT shapes. | `scenarios.npz`, `scenarios.json` |
| `05_matrix_imaging` | For each scenario: phase 1 picks where the BDT high is, phase 2 images it; static designs for comparison. Scores location error and coverage of the BDT top. Noise-free. | `results_05_mi.csv` |


