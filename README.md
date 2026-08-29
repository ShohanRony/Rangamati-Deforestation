# Rangamati Deforestation Monitoring

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22108374.svg)](https://doi.org/10.5281/zenodo.22108374)

Open-source reproduction and extension of a B.Sc. thesis on deforestation
monitoring in Rangamati, Chittagong Hill Tracts, Bangladesh, using
Google Earth Engine and Python (Google Colab). Covers a 30-year Landsat
time series (1993-2023) across six epochs, with multi-algorithm
land-cover classification, spatially independent accuracy validation,
and transition analysis.

## Study area

The whole of Rangamati district (all ten upazilas — not a sub-district
subset), Chittagong Hill Tracts, Bangladesh. The Earth Engine boundary
asset used throughout the scripts carries the attribute label
`NAME_2 == 'Parbattya Chattagram'`; this is a naming artifact on the
source asset, not a scope change — the underlying polygon is the
Rangamati district boundary, and analysis is performed and reported at
the district level only (no upazila-level breakdown).

## Data

- Landsat 5 (1993, 1998, 2008), Landsat 7 (2003), and Landsat 8 (2018, 2023)
  Collection 2 Level 2 surface reflectance, dry-season (January-April)
  composites.
- ESA WorldCover 2021 as training-label reference.
- A 200-point manually interpreted validation set (Script 06a), verified
  against high-resolution imagery (Script 06b re-extracts classifier
  predictions at these points once interpretation is complete).

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
> Monitoring in Rangamati, Bangladesh (1993–2023)* [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.22108374

## Method

Five land-cover classes (Dense Forest, Degraded Forest/Jhum, Water,
Agriculture/Settlement, Bare Land) are mapped for each epoch using a
Random Forest classifier, with a single spatial-block holdout split (the
study area is split into 0.1-degree geographic blocks, each assigned
wholly to train or holdout via one seeded random draw per block, so
nearby pixels never leak across the split — this is a single holdout,
not k-fold cross-validation). Accuracy is assessed with confusion
matrices and Olofsson et al. (2014)
area-adjusted estimation. Multi-decadal change is characterized with
land-cover transition matrices between consecutive epochs and over the
full 1993-2023 period.

## Repository structure

```
gee-scripts/        Google Earth Engine JavaScript scripts (run in the
                     GEE Code Editor, in numeric order)
python-analysis/     Olofsson area-adjusted accuracy computation
notebooks/           Exploratory analysis and figure generation
data/                Validation point sets (CSV / GeoJSON / KML)
```

### GEE scripts

| Script | Purpose |
|---|---|
| `01_data_exploration.js` | NDVI comparison between Landsat 5 and Landsat 8 epochs |
| `01_study_area_exploration.js` | Study area boundary and basic data exploration |
| `02_data_pipeline.js` | Multi-year Landsat composite pipeline with spectral indices |
| `03_classification.js` | Training sample collection via stratified sampling |
| `04_multi-algorithm_land-cover_classification.js` | Random Forest / SVM / CART classifier comparison with spatial-block validation |
| `05_Export_Classified_Maps.js` | Export classified maps, forest masks, and composites to Drive |
| `06a_Independent_visual_validation_points.js` | Generate the 200-point independent validation set |
| `06b_Re-classify_independent_validation.js` | Extract classifier predictions at validation points |
| `07_class_area_calculation.js` | Per-class area calculation for Olofsson accuracy |
| `08_sensitivity_analysis.js` | Forest-definition sensitivity analysis across NDVI thresholds |
| `09_Transition_Analysis.js` | Land-cover transition matrices between epochs and cumulative 1993-2023 |
| `10_data_quality_check.js` | Per-epoch Landsat scene count / cloud cover / valid-pixel diagnostic |
| `11_Spectral_Separability_Diagnostic.js` | Dense Forest vs Degraded/Jhum NDVI class-separability (M-statistic) per epoch, testing the class-boundary-instability hypothesis for the 2008 dip |

Each script is self-contained and can be run directly in the GEE Code
Editor; later scripts reuse the same preprocessing and classifier-build
logic as earlier ones (noted in each script's header).

### Python analysis

`python-analysis/olofsson_accuracy.py` computes area-adjusted accuracy
estimates and confidence intervals (Olofsson et al., 2014) from the
confusion matrices and class-area totals exported by the GEE scripts.

`notebooks/publication_figures.ipynb` (run in Google Colab) generates
every publication figure — study-area map, per-epoch classifications,
the 1993-vs-2023 comparison, the cumulative forest-change map, the NDVI
change map, and the forest-area trend chart — from the GeoTIFFs and CSVs
exported by the GEE scripts. No GIS desktop software is used anywhere
in this pipeline.

## Researcher

Shohinur Pervez Shohan, B.Sc. Computer Science & Engineering, Rangamati
Science and Technology University. Supervisor: Dhiman Sarma, Associate
Professor, RMSTU.
