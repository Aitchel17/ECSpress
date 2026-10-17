# ECSpress

Two-photon analysis of cerebral vessels and the perivascular space: vessel and PVS
width over arousal states, the polar partition of the PVS boundary, Radon anisotropy,
and PIV of the tissue around the vessel.

## Entry scripts, in running order

The scripts at the root are the entry points. The letter says which layer a script
reads; scripts under one letter are independent of each other.

| | reads | script | writes |
|---|---|---|---|
| a | one session folder | `a_analysis_pvs` | `paxfwhm.mat`, `polarcluster.mat`, `radon_result.mat`, `roilist.mat` under the session |
| a | one session folder | `a_check_pivquality` | `pivquality_pass.mat`, one mouse's sessions judged for PIV by eye |
| a | one session folder | `a_analysis_piv` | `pivtransition.mat`, `pivpolar.mat`; needs `a_analysis_pvs` and `a_check_pivquality` |
| b | one session's products | `b_reviewsession` | the review figures, and the grade through `ledger_quality` |
| c | the whole data tree | `c_dirtable` | each cohort's directory workbook and the joined dataset sheet |
| d | the sheet and every session's products | `d_centralize` | `centralized_<product>.mat`, then `centralized_bout.mat` and the series split by bout |
| e | the centralized state files | `06_analysis_integration/tablegeneration_main`, `_drift`, `_fwhmrelation` (scripts for now; they become functions) | the tables the figure scripts under `07_figuregeneration` read |

Every script opens by adding `09_dirstruct` and calling `dirs_ecspath`, which puts this
tree and its two siblings on the path. `archive_session` is the scheduled copy to the
archive disk and is run by name from a Windows task.

## Layout

| folder | what |
|---|---|
| `00_plotting` | figure helpers shared by the analysis classes and the figure scripts |
| `01_loaddata` | `ECSSession`, the session folder as an object |
| `02_preprocess` | stack downsampling and filtering |
| `03_roi` | ROI drawing and storage |
| `04_analysis` | the analysis classes: FWHM, polar cluster, Radon, PIV, event detection |
| `05_stateanalysis` | bout tables from the scorings, and `state_1d` |
| `06_analysis_integration` | the centralize and bout functions `d_centralize` calls, and the table helpers |
| `07_figuregeneration` | the figure scripts |
| `09_dirstruct` | the data-tree scan, the directory tables, the path and root resolution |
| `10_ledger` | the session quality register |
| `functions` | small helpers called from more than one place |
| `zz_notinuse` | retired code, kept off the path |
