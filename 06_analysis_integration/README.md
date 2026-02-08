# Secondary Analysis Repository
**Path**: `g:\03_program\03_secondary_analysis`

This repository handles the aggregation and secondary analysis of data processed by `01_ecspress` and `02_othersignal`.

## Dependencies
This repository **depends** on the following:
1. `01_ecspress`: For core image processing functions (`xy2heatmap`, etc.) and classes.
2. `02_othersignal`: For physiological signal analysis.

## Setup
Before running any script, you **MUST** run:
```matlab
setup_workspace
```
This script adds the necessary dependency folders to your MATLAB path.

## Structure
- `scripts/`: Analysis logic (moved from `ecspress/second_analysis`).
- `data/`: Local storage for aggregated `.mat` files.
