// ============================================================================
// Script 11 - Dense Forest vs Degraded/Jhum Spectral Separability Diagnostic
// ============================================================================
// Purpose: test the second (untested) hypothesis behind the 2003-2008-2018
// swing reported in Results Section 3.4 - that 2008 shows *worse spectral
// separability* between Dense Forest and Degraded/Jhum than other epochs,
// which would support "class-boundary instability" as a real contributing
// cause (as opposed to a confirmed deforestation/regrowth event).
//
// Script 10 already ruled out image quality (cloud cover, scene count,
// valid-pixel fraction) as the cause. This script complements that by
// directly measuring class separability per epoch using NDVI.
//
// Method: for each epoch, sample NDVI from pixels classified as
// Dense Forest (class 1) and Degraded/Jhum (class 2) in that epoch's
// RF-classified map, then compute the M-statistic (Kaufman & Remer, 1994),
// a standard remote-sensing class-separability measure:
//
//     M = |mean_1 - mean_2| / (stdDev_1 + stdDev_2)
//
// M > 1   -> well separated (histograms barely overlap)
// M < 1   -> poor separation (substantial overlap) -> supports the
//            class-boundary-instability hypothesis for that epoch
//
// ---------------------------------------------------------------------------
// SETUP
// ---------------------------------------------------------------------------
// studyArea and the NDVI-composite builder below are ported verbatim from
// gee-scripts/02_data_pipeline.js (the live script in this account) so this
// diagnostic uses exactly the same per-epoch composite as the rest of the
// pipeline: same cloud filter (CLOUD_COVER < 50, scene-level, no per-pixel
// QA_PIXEL masking - Script 02 itself doesn't apply that either), same
// Jan-Apr compositing window, same Roy et al. (2016) cross-sensor
// harmonization for Landsat 8. classifiedAssets below point at the 6
// RF_Classified_<year> images uploaded as EE assets under this project.

var studyArea = ee.FeatureCollection('projects/crypto-hallway-405211/assets/BGD_adm2')
  .filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'))
  .geometry();

var epochs = [1993, 1998, 2003, 2008, 2018, 2023];

var classifiedAssets = {
  1993: 'projects/crypto-hallway-405211/assets/RF_Classified_1993',
  1998: 'projects/crypto-hallway-405211/assets/RF_Classified_1998',
  2003: 'projects/crypto-hallway-405211/assets/RF_Classified_2003',
  2008: 'projects/crypto-hallway-405211/assets/RF_Classified_2008',
  2018: 'projects/crypto-hallway-405211/assets/RF_Classified_2018',
  2023: 'projects/crypto-hallway-405211/assets/RF_Classified_2023',
};

var DENSE_FOREST_CLASS = 1;
var DEGRADED_JHUM_CLASS = 2;
var SAMPLES_PER_CLASS = 500;   // per epoch, per class
var SEED = 42;

// ---------------------------------------------------------------------------
// NDVI composite builder - ported verbatim from 02_data_pipeline.js
// ---------------------------------------------------------------------------

function scaleL5(img) {
  var optical = img.select('SR_B.*').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true).copyProperties(img, img.propertyNames());
}
function scaleL7(img) {
  var optical = img.select('SR_B.*').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true).copyProperties(img, img.propertyNames());
}
function scaleL8(img) {
  var optical = img.select('SR_B.*').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true).copyProperties(img, img.propertyNames());
}

function harmonizeL8toL5(img) {
  var slopes = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]);
  var intercepts = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]);
  var harmonized = img
    .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'])
    .multiply(slopes).add(intercepts)
    .rename(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']);
  return img.addBands(harmonized, null, true).copyProperties(img, img.propertyNames());
}

// L5/L7 bands: B1=Blue, B2=Green, B3=Red, B4=NIR, B5=SWIR1, B7=SWIR2
function addNdviL5(img) {
  return img.addBands(img.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI'));
}
function addNdviL8(img) {
  // After harmonization, L8 bands are already renamed to L5-style naming.
  return img.addBands(img.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI'));
}

function makeL5Composite(year) {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL5).map(addNdviL5)
    .median().clip(studyArea).set('year', year);
}
function makeL7Composite(year) {
  return ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL7).map(addNdviL5)  // L7 same band structure as L5
    .median().clip(studyArea).set('year', year);
}
function makeL8Composite(year) {
  var start = (year === 2013) ? year + '-04-01' : year + '-01-01';
  var end   = (year === 2013) ? year + '-07-31' : year + '-04-30';
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL8).map(harmonizeL8toL5).map(addNdviL8)
    .median().clip(studyArea).set('year', year);
}

// Sensor-per-epoch mapping, exactly as Script 02:
// L5: 1993, 1998, 2008 | L7: 2003 | L8: 2018, 2023 (2013 excluded from study)
function buildNdviComposite(year) {
  if (year === 1993 || year === 1998 || year === 2008) return makeL5Composite(year);
  if (year === 2003) return makeL7Composite(year);
  if (year === 2018 || year === 2023) return makeL8Composite(year);
  throw new Error('No composite builder for year ' + year);
}

// ---------------------------------------------------------------------------
// Core routine
// ---------------------------------------------------------------------------

function computeSeparability(year) {
  var classified = ee.Image(classifiedAssets[year]);
  var ndvi = buildNdviComposite(year).select('NDVI');

  var denseMask = classified.eq(DENSE_FOREST_CLASS);
  var degradedMask = classified.eq(DEGRADED_JHUM_CLASS);

  var denseNdvi = ndvi.updateMask(denseMask).rename('ndvi_dense');
  var degradedNdvi = ndvi.updateMask(degradedMask).rename('ndvi_degraded');

  var denseSample = denseNdvi.sample({
    region: studyArea,
    scale: 30,
    numPixels: SAMPLES_PER_CLASS,
    seed: SEED,
    geometries: false,
  });

  var degradedSample = degradedNdvi.sample({
    region: studyArea,
    scale: 30,
    numPixels: SAMPLES_PER_CLASS,
    seed: SEED,
    geometries: false,
  });

  var denseStats = denseSample.reduceColumns({
    reducer: ee.Reducer.mean().combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true }),
    selectors: ['ndvi_dense'],
  });

  var degradedStats = degradedSample.reduceColumns({
    reducer: ee.Reducer.mean().combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true }),
    selectors: ['ndvi_degraded'],
  });

  var meanDense = ee.Number(denseStats.get('mean'));
  var stdDense = ee.Number(denseStats.get('stdDev'));
  var meanDegraded = ee.Number(degradedStats.get('mean'));
  var stdDegraded = ee.Number(degradedStats.get('stdDev'));

  var mStatistic = meanDense.subtract(meanDegraded).abs()
    .divide(stdDense.add(stdDegraded));

  return ee.Feature(null, {
    epoch: year,
    n_dense: denseSample.size(),
    n_degraded: degradedSample.size(),
    mean_ndvi_dense: meanDense,
    std_ndvi_dense: stdDense,
    mean_ndvi_degraded: meanDegraded,
    std_ndvi_degraded: stdDegraded,
    m_statistic: mStatistic,
  });
}

var results = ee.FeatureCollection(epochs.map(computeSeparability));

// Console preview
print('Spectral separability (M-statistic) by epoch', results);
print('Reminder: M < 1 = poor separation (overlap); M > 1 = well separated. ' +
      'If 2008 is the lowest M-statistic among all six epochs, that supports ' +
      'class-boundary instability as a contributing cause of the 2003-2008-2018 swing.');

// Export for the paper (Section 3.4 supplementary evidence / Discussion)
Export.table.toDrive({
  collection: results,
  description: 'DenseForest_Degraded_Separability_AllEpochs',
  fileFormat: 'CSV',
  folder: 'Rangamati_Deforestation',
});

// ---------------------------------------------------------------------------
// RESULTS (run 2026-08-28, against the corrected v2 classification pipeline)
// ---------------------------------------------------------------------------
// epoch | m_statistic | mean_ndvi_dense | mean_ndvi_degraded | n_dense | n_degraded
// 1993  | 1.173       | 0.660           | 0.514               | 334     | 71
// 1998  | 0.757       | 0.620           | 0.485               | 346     | 66
// 2003  | 0.999       | 0.640           | 0.478               | 347     | 71
// 2008  | 0.921       | 0.591           | 0.480               | 171     | 195
// 2018  | 1.270       | 0.698           | 0.564               | -       | 71
// 2023  | 1.148       | 0.688           | 0.558               | -       | -
//
// 2008 is NOT the single lowest M-statistic (1998 is lower, at 0.757), so this
// does not cleanly confirm 2008-specific class-boundary instability. Three
// epochs (1998, 2003, 2008) all sit at or below M~1 while 1993/2018/2023 are
// clearly well-separated (M>1.14). Notably, 2008's Dense/Degraded sample split
// (171/195) is a near-inversion of every other epoch's ~330/70 split, which is
// a distinct finding from the spectral-separability question itself.
