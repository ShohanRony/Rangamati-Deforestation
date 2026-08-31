# Rangamati Land-Cover Change, 1993-2023

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22108374.svg)](https://doi.org/10.5281/zenodo.22108374)

Open-source reproduction and extension of a B.Sc. thesis on land-cover
monitoring in Rangamati District, Chittagong Hill Tracts, Bangladesh,
using Google Earth Engine and Python (Google Colab). The pipeline
reconstructs land cover across six Landsat epochs (1993-2023), compares
three classifiers under a spatially blocked evaluation design, assesses
the 2023 map with area-adjusted accuracy estimation, tests the
sensitivity of forest-area estimates to a single operational threshold,
and maps land-cover transitions between epochs.

This is a single-district case study, not a survey of the full
Chittagong Hill Tracts (three districts); see "Study area" below.

## Study area and objectives

Rangamati District, the largest of the three Chittagong Hill Tracts
districts, all ten upazilas, no upazila-level breakdown. The Earth Engine
boundary asset carries the attribute label `NAME_2 == 'Parbattya
Chattagram'` — a naming artifact on the source asset, not a scope change;
the underlying polygon is the Rangamati district boundary, and every
result in this repository is reported at the district level only.

Objectives:

1. Reconstruct five-class land cover for six Landsat epochs under one
   documented processing framework, including inter-sensor harmonisation.
2. Compare Random Forest, SVM, and CART under a spatial-block holdout
   rather than a random split, to avoid inflating accuracy through
   spatially autocorrelated train/test leakage.
3. Assess the 2023 map with an independent, manually interpreted
   reference sample and Olofsson et al. (2014) area-adjusted estimation.
4. Test how sensitive the reported forest area is to the NDVI threshold
   that separates Dense Forest from Degraded Forest/Jhum.
5. Map land-cover transitions between consecutive epochs and over the
   full 1993-2023 period, and report — rather than resolve — an
   unexplained 2003-2008-2018 forest-composition swing.

### Boundary data source

The Bangladesh Level-2 administrative boundary (`BGD_adm2`) used to
derive the study-area polygon is a [GADM](https://gadm.org) dataset,
mirrored for download by [DIVA-GIS](https://www.diva-gis.org/). It is
**not redistributed in this repository** — GADM's license permits free
academic/non-commercial use but not redistribution without prior
permission. To reproduce the boundary asset used in the scripts:

1. Download the Bangladesh Level-2 administrative boundary shapefile
   from [DIVA-GIS](https://www.diva-gis.org/data.html) or directly from
   [GADM](https://gadm.org/download_country.html) (select Bangladesh).
2. Upload the shapefile to Google Earth Engine as a table asset named
   `BGD_adm2`.
3. In the scripts, the study-area polygon is selected with
   `NAME_2 == 'Parbattya Chattagram'` (see "Study area" above).

## Data

- **Six Landsat epochs**: Landsat 5 (1993, 1998, 2008), Landsat 7 (2003),
  and Landsat 8 (2018, 2023), USGS Collection 2 Level 2 surface
  reflectance, dry-season (1 January - 30 April) composites. A 2013
  Landsat 8 epoch was explored but excluded from the production analysis
  (insufficient cloud-free coverage in the compositing window).
- **Five operational land-cover classes**: Dense Forest, Degraded
  Forest/Jhum, Water, Agriculture/Settlement, Bare Land — derived from
  ESA WorldCover 2021, with tree cover split into Dense Forest
  (NDVI >= 0.55) and Degraded Forest/Jhum (0.25 <= NDVI < 0.55) plus
  WorldCover shrubland/grassland. These are operational remote-sensing
  categories defined by a spectral rule, not field measurements of
  canopy condition (see "Limitations" below).
- ESA WorldCover 2021 as the training-label reference.
- A 200-point manually interpreted validation set (Script 06a), verified
  against high-resolution imagery (Script 06b re-extracts classifier
  predictions at these points once interpretation is complete).

## Workflow

All processing is implemented in the Google Earth Engine JavaScript API,
with area-adjusted accuracy computed separately in Python and every
publication figure generated in a Colab notebook — no desktop GIS
software is used anywhere in the pipeline.

### Preprocessing and sensor harmonisation

Each Landsat scene is masked at pixel level (`QA_PIXEL`: dilated cloud,
cloud, cloud shadow, snow; `QA_RADSAT`: saturation) and scaled to surface
reflectance (`DN * 0.0000275 - 0.2`), then combined into a per-epoch
median composite over the dry-season window. Landsat 8 reflectance is
harmonised to the Landsat 5/7 spectral domain using the OLI-to-ETM+
continuity coefficients of Roy et al. (2016) before spectral indices
(NDVI, NDWI, EVI, NBR) and terrain predictors (SRTM elevation and slope)
are computed, giving a twelve-band feature stack per epoch.

### Classifier comparison and spatial-block holdout

Random Forest (500 trees), SVM (RBF kernel), and CART are compared under
a single **spatial-block holdout**, not k-fold cross-validation: the
study area is divided into 0.1-degree geographic blocks, and each block
is assigned wholly to the training set or the holdout set via one seeded
random draw per block (seed 42) — so nearby pixels never leak across the
split, which would otherwise inflate reported accuracy. Random Forest is
carried forward into production mapping.

### 2023 area-adjusted accuracy assessment

The 2023 map is assessed against 200 manually interpreted reference
points using the design-based estimator of Olofsson et al. (2014), which
converts the raw confusion matrix into area-adjusted overall accuracy,
per-class user's/producer's accuracies, and confidence intervals,
weighted by each class's mapped area. This assessment applies **only**
to the 2023 map (see "Limitations").

### NDVI-threshold sensitivity analysis

Because the Dense Forest / Degraded Forest/Jhum boundary is a single
NDVI threshold (operationally fixed at 0.55), the 2023 classification is
re-run at 0.45 and 0.60 to test how much the reported forest area
depends on this one parameter, independent of classifier choice.

### M-statistic diagnostic

A supplementary diagnostic tests whether the six epochs differ in how
well Dense Forest and Degraded Forest/Jhum separate spectrally, using
the M-statistic of Kaufman & Remer (1994):
`M = |mean_1 - mean_2| / (stdDev_1 + stdDev_2)` computed per epoch from
NDVI sampled at classified pixels. This diagnostic intentionally reuses
an earlier compositing routine (Script 02) rather than the production
pipeline (Script 04/05) and is reported as exploratory supplementary
evidence, not a formal accuracy assessment (see Script 11's header).

## Principal verified results

Figures below are computed by the scripts in this repository and
reported in the accompanying manuscript. Historical-epoch figures are
unadjusted, mapped reconstructions (see "Limitations"); only the 2023
figures are area-adjusted.

| Result | Value |
|---|---|
| Classifier comparison, spatial-block holdout | Random Forest 77.10%, SVM 78.04%, CART 72.90% overall accuracy |
| 2023 area-adjusted overall accuracy | 73.17% (95% CI 64.34-82.00%) |
| Combined forest area change, 1993 -> 2023 | +20,693 ha (+4.4%) |
| Dense Forest area change, 1993 -> 2023 | -8,235 ha |
| NDVI-threshold sensitivity (0.45 vs 0.60) | Dense Forest area changes ~27%; combined forest area changes ~0.7% |
| 2003-2008-2018 forest-composition swing | Not explained by input-image quality (Script 10); M-statistic diagnostic (Script 11) neither confirms nor excludes class-boundary instability; reported as unresolved |

## Limitations and interpretation boundaries

- **Area-adjusted validation applies only to the 2023 map.** Historical
  epochs (1993, 1998, 2003, 2008, 2018) are unadjusted, mapped-area
  reconstructions produced by temporal transfer of the 2023-trained
  classifier, not by independent per-epoch retraining or validation.
- **The 200-point reference sample is not fully independent** of the
  map: point locations were selected during an earlier classification
  iteration and retained after a later classifier correction (see the
  caveat in `scripts/python/olofsson_accuracy.py`).
- **Dense Forest and Degraded Forest/Jhum are operational remote-sensing
  categories**, defined by a single NDVI threshold against WorldCover
  tree cover — not direct field measurements of forest condition.
- **The 2003-2008-2018 swing is reported as unresolved**, not as
  confirmed deforestation, degradation, or recovery. Two diagnostics
  (Scripts 10 and 11) narrow but do not identify the cause.
- No causal or policy claim is made anywhere in this pipeline; outputs
  describe mapped/estimated land-cover composition and change only.

Full detail, including which script produces which caveat, is in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md), section 5.

## Data availability

The processed outputs of this study — classified land-cover GeoTIFFs
(1993-2023), forest masks, the NDVI change raster, the true-color
composite, validation points, and the publication figures — are archived
on Zenodo:

| | |
|---|---|
| **DOI** | [10.5281/zenodo.22108374](https://doi.org/10.5281/zenodo.22108374) |
| **Resource type** | Dataset |
| **Publisher** | Zenodo |
| **Language** | English |
| **License** | CC BY 4.0 |

### Citation

> Shohan, S. P., & Sarma, D. (2026). *Landsat-Based Random Forest Land
> Cover Classification and Validation Dataset for Deforestation
> Monitoring in Rangamati, Bangladesh (1993-2023)* [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.22108374

## License

**No license has been selected for the code in this repository yet** —
see [`LICENSE`](LICENSE) for a placeholder and common options. The
CC BY 4.0 license above covers the Zenodo *data* archive only, not the
source code here.

## Repository structure

```
scripts/gee/         Google Earth Engine JavaScript scripts (run in the
                      GEE Code Editor, in numeric order — see
                      docs/REPRODUCIBILITY.md for the exact sequence)
scripts/python/       Olofsson area-adjusted accuracy computation
notebooks/            Exploratory analysis and figure generation
data/                 Validation point sets (CSV / GeoJSON / KML) and the
                      Script 11 spectral-separability diagnostic output
docs/                 Reproducibility guide (assets, execution order,
                      dependencies, validated-vs-unadjusted results)
LICENSE               Placeholder — no license chosen yet (see above)
```

### GEE scripts

Each script's header states its status (production vs. exploratory),
purpose, inputs, outputs, and dependencies in full; this table
summarises them. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)
for the exact execution order and required assets.

| Script | Status | Purpose |
|---|---|---|
| `01_data_exploration.js` | exploratory | Earliest NDVI comparison, Landsat 5 vs 8, bounding-rectangle approximation |
| `01_study_area_exploration.js` | exploratory | Real boundary preview + true-color composites |
| `02_data_pipeline.js` | exploratory | Early full-pipeline draft; different harmonisation coefficients, includes 2013, naive NDVI differencing |
| `03_classification.js` | exploratory | Training-sample prototype (superseded by Script 04's inline sampling) |
| `04_multi-algorithm_land-cover_classification.js` | **production** | Random Forest / SVM / CART comparison under spatial-block holdout |
| `05_Export_Classified_Maps.js` | **production** | Export classified maps, forest masks, and composites to Drive |
| `06a_Independent_visual_validation_points.js` | **production** | Generate the 200-point independent validation set |
| `06b_Re-classify_independent_validation.js` | **production** | Extract classifier predictions at validation points |
| `07_class_area_calculation.js` | **production** | Per-class area calculation for Olofsson accuracy |
| `08_sensitivity_analysis.js` | **production** | Forest-definition sensitivity analysis across NDVI thresholds |
| `09_Transition_Analysis.js` | **production** | Land-cover transition matrices between epochs and cumulative 1993-2023 |
| `10_data_quality_check.js` | **production** | Per-epoch Landsat scene count / cloud cover / valid-pixel diagnostic |
| `11_Spectral_Separability_Diagnostic.js` | **production** (supplementary) | Dense Forest vs Degraded/Jhum NDVI class-separability (M-statistic) per epoch |

### Python analysis

`scripts/python/olofsson_accuracy.py` computes area-adjusted accuracy
estimates and confidence intervals (Olofsson et al., 2014) from the
confusion matrix and class-area totals exported by the GEE scripts.

### Notebook

`notebooks/publication_figures.ipynb` (run in Google Colab) generates
every publication figure — study-area map, per-epoch classifications,
the 1993-vs-2023 comparison, the cumulative forest-change map, the NDVI
change map, and the forest-area trend chart — from the GeoTIFFs and CSVs
exported by the GEE scripts.

## Reproducibility

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for:

- required Google Earth Engine assets and how to create them under your
  own Cloud project;
- exact script execution order, including which steps are independent
  of each other;
- Python dependencies;
- expected inputs and outputs for every stage;
- which results are area-adjusted/validated and which are unadjusted
  mapped reconstructions;
- known limitations and private-data restrictions.

## Researcher

Shohinur Pervez Shohan, B.Sc. Computer Science & Engineering, Rangamati
Science and Technology University. Supervisor: Dhiman Sarma, Associate
Professor, RMSTU.
