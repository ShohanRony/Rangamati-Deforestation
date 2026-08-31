# Reproducibility guide

This document lists, in execution order, everything needed to reproduce
every numerical result reported for the Rangamati District study. It
complements the top-level README; script-by-script purpose/input/output
detail lives in each script's own header comment.

## 1. Required Google Earth Engine assets

All production scripts (04 onward) reference **private** assets under one
GEE Cloud project. To reproduce them under your own project:

| Asset (as used in scripts) | Source | How to create it |
|---|---|---|
| `.../assets/BGD_adm2` | GADM Bangladesh Level-2 boundary, mirrored by [DIVA-GIS](https://www.diva-gis.org/data.html) | Download shapefile, upload to GEE as a Table asset named `BGD_adm2`. See README "Boundary data source" for the license note. |
| `.../assets/ValidationPoints_200` (or `_Reupload`) | Produced by Script 06a, manually interpreted, then uploaded from `data/ValidationPoints_200_ForReupload.geojson` | Upload the geojson in `data/` as a Table asset; see Script 06b header for the exact steps. |
| `.../assets/RF_Classified_<year>` (x6, 1993/1998/2003/2008/2018/2023) | Exported by Script 05 as GeoTIFF, then re-uploaded as Image assets | Run Script 05 first, download the GeoTIFFs from Drive, upload each as an Image asset named `RF_Classified_<year>`. Only required for Script 11 (spectral-separability diagnostic); not required for Scripts 04-09. |

Every script's `projects/crypto-hallway-405211/assets/...` path is specific
to the original author's GEE Cloud project. To run any script, replace the
project segment with your own project ID after creating the assets above.

## 2. Execution order

```
scripts/gee/03_classification.js                 (optional; exploratory prototype, not required downstream)
scripts/gee/04_multi-algorithm_land-cover_classification.js   -> classifier comparison (Table: OA/kappa per classifier)
scripts/gee/05_Export_Classified_Maps.js          -> 6x classified GeoTIFF, forest masks, composites, NDVI-change stack (Drive)
scripts/gee/06a_Independent_visual_validation_points.js -> ~200 candidate points for manual interpretation
   [manual step: interpret each point against high-resolution imagery -> manual_class]
scripts/gee/06b_Re-classify_independent_validation.js   -> confusion-matrix CSV (needs 06a's uploaded asset)
scripts/gee/07_class_area_calculation.js          -> per-class mapped-area CSV (2023)
scripts/python/olofsson_accuracy.py               -> area-adjusted OA, per-class UA/PA, CIs (needs 06b + 07 CSVs)
scripts/gee/08_sensitivity_analysis.js            -> NDVI-threshold sensitivity (independent of the above)
scripts/gee/09_Transition_Analysis.js             -> transition matrices, net/annualised change (independent of the above)
scripts/gee/10_data_quality_check.js              -> per-epoch scene/cloud/valid-pixel diagnostic (independent)
scripts/gee/11_Spectral_Separability_Diagnostic.js -> M-statistic per epoch (needs Script 05's 6 classified-map assets)
notebooks/publication_figures.ipynb               -> all publication figures (needs Script 05's Drive exports)
```

`scripts/gee/01_data_exploration.js`, `01_study_area_exploration.js`, and
`02_data_pipeline.js` are earlier exploratory scripts, not required for any
of the above and not part of the production pipeline (see their headers).

Each GEE script other than 06b, 07->olofsson, and 11 is **self-contained**:
it repeats its own preprocessing, sampling, and classifier training rather
than importing a previous script's in-session state. This is deliberate
(GEE Code Editor sessions do not persist between scripts) but means that
each script's stratified sample is regenerated (same seed, same result)
rather than reused — read literally, "training data" is rebuilt, not
carried over, at every script boundary except where a file/asset is
explicitly exported and re-imported (noted above).

## 3. Python dependencies

```
pandas
numpy
```

No `requirements.txt` is pinned to specific versions; `olofsson_accuracy.py`
uses only basic `pandas`/`numpy` operations (crosstab, arithmetic), so any
reasonably recent version of both should work. `publication_figures.ipynb`
is run in Google Colab, which provides its own pinned environment.

## 4. Expected inputs / outputs per stage

| Stage | Reads | Writes |
|---|---|---|
| Script 04 | GEE public collections + `BGD_adm2` | console only (OA/kappa) |
| Script 05 | same | Drive: 6x classified GeoTIFF, 6x forest mask, 2023 RGB composite, NDVI-change stack |
| Script 06a | same + `ESA/WorldCover/v200/2021` | FeatureCollection of candidate points (export as needed) |
| Script 06b | `ValidationPoints_200_*` asset (from 06a, manually labelled) | CSV: manual_class + rf_class per point |
| Script 07 | same as 04 | console (per-class ha) + `RF_MappedArea_PerClass_2023.csv` |
| `olofsson_accuracy.py` | `RF_MappedArea_PerClass_2023.csv`, `ValidationPoints_200_Reclassified_v2.csv` | console: OA, per-class UA/PA, CIs |
| Script 08 | same as 04 | console (Dense Forest ha per NDVI threshold) |
| Script 09 | same as 04 | console (transition matrices, ha + ha/yr) |
| Script 10 | Landsat collections only | console (scene count, cloud %, valid-pixel %) |
| Script 11 | 6x `RF_Classified_<year>` assets (from 05) | Drive CSV: M-statistic per epoch (`data/DenseForest_Degraded_Separability_AllEpochs.csv`) |

## 5. Validated vs. unadjusted results — read before citing a number

- **Area-adjusted validation applies only to the 2023 map.** The overall
  accuracy of 73.17% (95% CI 64.34-82.00%) and the per-class user's/
  producer's accuracies come from the Olofsson estimator applied to the
  200-point reference sample against the 2023 classification, computed by
  `olofsson_accuracy.py`.
- **The 200-point reference sample is not fully independent of the map**:
  point locations were selected during an earlier classification
  iteration and retained after a later classifier correction; see the
  caveat paragraph at the top of `olofsson_accuracy.py` and the
  manuscript's accuracy-assessment section for the full disclosure.
- **Historical epochs (1993, 1998, 2003, 2008, 2018) are unadjusted,
  mapped-area reconstructions**, not independently validated. They are
  produced by applying one classifier trained on the 2023 composite to
  each historical composite (temporal transfer), not by retraining per
  epoch. No equivalent 200-point (or any) manual reference set exists for
  these epochs.
- **Dense Forest / Degraded Forest/Jhum are operational NDVI-threshold
  classes**, not field measurements of canopy condition (see Script 04's
  header for the exact class definitions).
- **Script 11's M-statistic diagnostic is supplementary, not a formal
  accuracy assessment** — it deliberately reuses Script 02's (exploratory)
  compositing routine rather than the Script 04/05 production pipeline;
  see Script 11's own header for the exact deviation.
- **The 2003-2008-2018 forest-composition swing is reported as
  unresolved.** Script 10 rules out input-image quality as the cause;
  Script 11 does not confirm or exclude class-boundary instability. No
  script in this repository demonstrates a confirmed cause.

## 6. Known limitations and private-data restrictions

- The `BGD_adm2` boundary asset is derived from a GADM dataset that
  permits free academic/non-commercial use but not redistribution; it is
  not included in this repository (see README "Boundary data source").
- All GEE asset paths (`projects/crypto-hallway-405211/assets/...`) point
  to a private Cloud project. Reproducing the pipeline requires creating
  equivalent assets under your own project (Section 1, above).
- Classified GeoTIFFs, forest masks, the NDVI-change raster, and the
  publication figures are not stored in this repository; they are
  archived on Zenodo (see README "Data availability").
