# Notebooks

| Notebook | What it does | Main outputs |
|---|---|---|
| `01_velocity_model` | Builds the 3D Vp and Vs model from the known geology: layers, faults (Calibugan et al. 2022 structure) and greywacke basement. No BDT. | `velocity_model.npz` |
| `02_designs` | Station layouts and deployment schedule: Trillium rings, static node designs, and the phased design (wide reconnaissance grid, then dense rolling patches). Also the station bank phase 2 can be placed on. | `stations.csv`, `designs.json` |
| `03_travel_times` | P travel times from every station to every model cell (eikonal), plus 1D tables. Slow the first time. | `tt/` |
| `04_scenarios` | Possible "true" Earths: conductor positions sampled around the drill site, combined with several BDT shapes. | `scenarios.npz`, `scenarios.json` |
| `05_matrix_imaging` | For each scenario: phase 1 picks where the BDT high is, phase 2 images it; static designs for comparison. Scores location error and coverage of the BDT top. Noise-free. | `results_05_mi.csv` |
| `06_noise` | Repeats 05 with noise that shrinks with recording time, over a range of noise levels. | `results_06_noise.csv` |
| `07_model_error` | Repeats 05 with the real Earth differing from the model (velocities, clay thickness, greywacke depth). | `results_07_model_error.csv` |
| `08_hypocentres` | Earthquake location accuracy for the permanent network (GeoNet, plus Trilliums), for the future catalogue. | `results_08_hypo.csv` |
| `09_comparison` | Pulls 05 to 08 into one table and figure for choosing the design. | |
| `10_tomography` | Shows LET and attenuation tomography are too coarse to locate the BDT: ray coverage, checkerboards, BDT recovery. | |

`mi_common.py` holds the matrix imaging functions shared by 05, 06 and 07.


