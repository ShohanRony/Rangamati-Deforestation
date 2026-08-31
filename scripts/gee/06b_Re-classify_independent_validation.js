// ============================================================================
// Rangamati Land-Cover Change, 1993-2023
// Script 06b - Re-classify Independent Validation Points (PRODUCTION)
// ----------------------------------------------------------------------------
// STATUS: production. Output feeds scripts/python/olofsson_accuracy.py.
//
// Purpose:  The 200-point validation set (Script 06a) was manually
//           interpreted against 2023 high-resolution imagery, giving a
//           reference label (manual_class) per point. This script extracts
//           the corresponding RF-predicted class (rf_class) for each point
//           from the 2023 classified map, so the confusion matrix / Olofsson
//           area-adjusted accuracy can be computed. Manual reference labels
//           are produced separately in Script 06a and are not touched here.
// Inputs:   ValidationPoints_200_ForReupload.geojson (data/), uploaded as a
//           GEE Table asset (Assets tab > NEW > Table Upload); replace
//           VALIDATION_ASSET below with the resulting asset path, e.g.
//           'projects/crypto-hallway-405211/assets/ValidationPoints_200'.
// Outputs:  Table export (CSV) with manual_class and rf_class per point ->
//           scripts/python/olofsson_accuracy.py input.
// Depends:  Script 06a (validation-point geometries and manual labels).
// ============================================================================

var VALIDATION_ASSET = 'projects/crypto-hallway-405211/assets/ValidationPoints_200_Reupload';

// ---- STUDY AREA ----
var studyArea = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
).filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
var roi = studyArea.geometry();

// ---- PREPROCESSING (identical to Script 07) ----
function maskAndScale(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));
  var saturated = img.select('QA_RADSAT').eq(0);
  var optical = img.select('SR_B.*').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true)
    .updateMask(clear).updateMask(saturated)
    .copyProperties(img, img.propertyNames());
}

function harmonizeL8(img) {
  var slopes     = ee.Image.constant([0.8850,0.9317,0.9372,0.8339,0.8639,0.9165]);
  var intercepts = ee.Image.constant([0.0183,0.0123,0.0123,0.0448,0.0306,0.0116]);
  var l8 = img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']);
  return l8.multiply(slopes).add(intercepts)
    .rename(['Blue','Green','Red','NIR','SWIR1','SWIR2'])
    .copyProperties(img, img.propertyNames());
}

function addIndices(img) {
  var ndvi = img.normalizedDifference(['NIR','Red']).rename('NDVI');
  var ndwi = img.normalizedDifference(['Green','NIR']).rename('NDWI');
  var evi  = img.expression(
    '2.5 * (NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1)',
    {NIR: img.select('NIR'), RED: img.select('Red'), BLUE: img.select('Blue')}
  ).rename('EVI');
  var nbr  = img.normalizedDifference(['NIR','SWIR2']).rename('NBR');
  return img.addBands([ndvi, ndwi, evi, nbr]);
}

var srtm      = ee.Image('USGS/SRTMGL1_003').clip(roi);
var elevation = srtm.select('elevation');
var slope     = ee.Terrain.slope(srtm);
var featureBands = ['Blue','Green','Red','NIR','SWIR1','SWIR2',
                    'NDVI','NDWI','EVI','NBR','elevation','slope'];

// ---- BUILD 2023 COMPOSITE (identical to Script 07) ----
var composite2023 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2023-01-01', '2023-04-30')
  .map(maskAndScale)
  .map(harmonizeL8)
  .map(addIndices)
  .median()
  .clip(roi)
  .addBands(elevation.rename('elevation'))
  .addBands(slope.rename('slope'))
  .select(featureBands);

// ---- TRAINING DATA (identical to Script 07) ----
var worldcover = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(roi);
var ndvi2023   = composite2023.select('NDVI');

var dense    = worldcover.eq(10).and(ndvi2023.gte(0.55));
var degraded = worldcover.eq(10).and(ndvi2023.gte(0.25)).and(ndvi2023.lt(0.55))
               .or(worldcover.eq(20)).or(worldcover.eq(30));

var classMap = ee.Image(0)
  .where(dense,                                    1)
  .where(degraded,                                 2)
  .where(worldcover.eq(80),                        3)
  .where(worldcover.eq(40).or(worldcover.eq(50)),  4)
  .where(worldcover.eq(60),                        5)
  .rename('landcover')
  .updateMask(worldcover.gt(0));
classMap = classMap.updateMask(classMap.neq(0));

var samples = classMap.addBands(composite2023).stratifiedSample({
  numPoints: 150,
  classBand: 'landcover',
  region:    roi,
  scale:     30,
  seed:      42,
  geometries: true
});

// ---- Spatial block assignment (identical to Script 07) ----
samples = samples.map(function(f) {
  var c   = f.geometry().coordinates();
  var bx  = ee.Number(c.get(0)).multiply(10).floor();
  var by  = ee.Number(c.get(1)).multiply(10).floor();
  var bid = bx.format('%d').cat('_').cat(by.format('%d'));
  return f.set('block_id', bid);
});

var uniqueBlocks = samples.distinct('block_id').randomColumn('rand', 42);
var blockFolds = uniqueBlocks.map(function(b) {
  var isHoldout = ee.Number(b.get('rand')).gte(0.7);
  return b.set('block_fold', ee.Algorithms.If(isHoldout, 1, 0));
});

var joinFilter = ee.Filter.equals({leftField: 'block_id', rightField: 'block_id'});
var joined = ee.Join.saveFirst('blockMatch').apply(samples, blockFolds, joinFilter);
samples = joined.map(function(f) {
  var match = ee.Feature(f.get('blockMatch'));
  return f.set('block_fold', match.get('block_fold'));
});

var trainSet = samples.filter(ee.Filter.eq('block_fold', 0));

// ---- TRAIN RF (identical to Script 07 — this is the classifier that
//      produced the 396,367 / 99,640 / 30,624 / 35,084 / 14,820 ha areas) ----
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({
    features:        trainSet,
    classProperty:   'landcover',
    inputProperties: featureBands
  });

var classified2023 = composite2023.classify(rf).rename('rf_class_new');

// ================================================================
// RE-CLASSIFY THE 200 INDEPENDENT VALIDATION POINTS
// ================================================================
var validationPoints = ee.FeatureCollection(VALIDATION_ASSET);

var reclassified = classified2023.reduceRegions({
  collection: validationPoints,
  reducer: ee.Reducer.first(),
  scale: 30
});

// reduceRegions() with ee.Reducer.first() on a renamed single-band image
// still names its output column 'first', not the band name. Rename it
// explicitly so downstream code and the exported CSV have a clear,
// unambiguous column name.
reclassified = reclassified.map(function(f) {
  return f.set('rf_class_new', f.get('first'));
});

print('Re-classified validation points (first 5):', reclassified.limit(5));

// ---- Quick console confusion check (manual_class vs rf_class_new) ----
var confusion = reclassified.errorMatrix('manual_class', 'rf_class_new');
print('=== Independent validation — RAW confusion matrix ===');
print('Raw OA:', confusion.accuracy());
print('Kappa:', confusion.kappa());
print('Matrix:', confusion);

// ---- Export merged table for the Python Olofsson recomputation ----
Export.table.toDrive({
  collection: reclassified,
  description: 'ValidationPoints_200_Reclassified_v2',
  fileFormat: 'CSV',
  folder: 'Rangamati_Deforestation'
});

print('Export submitted: ValidationPoints_200_Reclassified_v2.csv');
print('This CSV has manual_class (unchanged reference) + rf_class_new (updated prediction).');
print('Use it, together with the new class areas from Script 07, to redo the Olofsson computation.');
