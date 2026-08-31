// ============================================================================
// Rangamati Land-Cover Change, 1993-2023
// Script 05 - Export Classified Maps (PRODUCTION)
// ----------------------------------------------------------------------------
// STATUS: production. Preprocessing/harmonisation/classifier setup mirrors
//   Script 04 exactly (same seed, same block split, same six epochs).
//
// Purpose:  Export the six per-epoch classified maps and supporting rasters
//           to Google Drive for downstream area calculation (Script 07),
//           transition analysis (Script 09), and figure generation
//           (notebooks/publication_figures.ipynb).
// Inputs:   projects/crypto-hallway-405211/assets/BGD_adm2 (private asset,
//           see docs/REPRODUCIBILITY.md); ESA/WorldCover/v200/2021;
//           Landsat 5/7/8 C02 T1_L2 collections.
// Outputs:  Drive exports -
//             1. RF-classified maps, GeoTIFF, one per epoch (6 files)
//             2. 2023 RGB composite (basemap reference)
//             3. NDVI change stack (visualization)
//             4. Forest binary mask, one per epoch (Python area calc.)
// Depends:  none (self-contained; regenerates its own training sample).
// Params:   seed=42; compositing window 1 Jan - 30 Apr (consistent with
//           Scripts 07/08/09).
// ============================================================================

// ---- STUDY AREA ----
var studyArea = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
).filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));

var roi = studyArea.geometry();

// ---- PREPROCESSING FUNCTIONS (same as Script 04) ----
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

function harmonizeL5(img) {
  var optical = img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']);
  return optical.rename(['Blue','Green','Red','NIR','SWIR1','SWIR2'])
    .copyProperties(img, img.propertyNames());
}

function harmonizeL7(img) {
  var optical = img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']);
  return optical.rename(['Blue','Green','Red','NIR','SWIR1','SWIR2'])
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

// SRTM terrain
var srtm      = ee.Image('USGS/SRTMGL1_003').clip(roi);
var elevation = srtm.select('elevation');
var slope     = ee.Terrain.slope(srtm);

// ---- BUILD COMPOSITE PER EPOCH ----
var epochs = [
  {year: 1993, sensor: 'L5', col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5},
  {year: 1998, sensor: 'L5', col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5},
  {year: 2003, sensor: 'L7', col: 'LANDSAT/LE07/C02/T1_L2', harmonize: harmonizeL7},
  {year: 2008, sensor: 'L5', col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5},
  {year: 2018, sensor: 'L8', col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8},
  {year: 2023, sensor: 'L8', col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8},
];

var featureBands = ['Blue','Green','Red','NIR','SWIR1','SWIR2',
                    'NDVI','NDWI','EVI','NBR','elevation','slope'];

function buildComposite(ep) {
  var col = ee.ImageCollection(ep.col)
    .filterBounds(roi)
    .filterDate(ep.year + '-01-01', ep.year + '-04-30')
    .map(maskAndScale)
    .map(ep.harmonize)
    .map(addIndices);
  var median = col.median().clip(roi);
  return median
    .addBands(elevation.rename('elevation'))
    .addBands(slope.rename('slope'))
    .select(featureBands);
}

// ---- REBUILD CLASSIFIER (RF) ----
// Training samples — same logic as Script 04
var worldcover = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(roi);
var ndvi2023   = buildComposite({
  year: 2023, sensor: 'L8',
  col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8
}).select('NDVI');

var treeMask     = worldcover.eq(10);
var dense        = treeMask.and(ndvi2023.gte(0.55));
var degraded     = treeMask.and(ndvi2023.gte(0.25)).and(ndvi2023.lt(0.55))
                   .or(worldcover.eq(20)).or(worldcover.eq(30));

var classMap = ee.Image(0)
  .where(dense,              1)
  .where(degraded,           2)
  .where(worldcover.eq(80),  3)
  .where(worldcover.eq(40).or(worldcover.eq(50)), 4)
  .where(worldcover.eq(60),  5)
  .rename('landcover')
  .updateMask(worldcover.gt(0));
classMap = classMap.updateMask(classMap.neq(0));

var composite2023 = buildComposite({
  year: 2023, sensor: 'L8',
  col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8
});

var samples = classMap.addBands(composite2023).stratifiedSample({
  numPoints: 150,
  classBand: 'landcover',
  region: roi,
  scale: 30,
  seed: 42,
  geometries: true
});

// ---------------------------------------------------------------------------
// Spatial block assignment (see Script 04 for full explanation)
// ---------------------------------------------------------------------------
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

print('Distinct blocks:', ee.FeatureCollection(uniqueBlocks).size());
print('Training samples (this export run):', trainSet.size());

// Train RF classifier
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({features: trainSet, classProperty: 'landcover', inputProperties: featureBands});

print('Classifier trained — starting exports...');

// ================================================================
// PART 1: EXPORT CLASSIFIED MAPS (6 epochs)
// ================================================================

// Color palette for visualization
var classVis = {
  min: 1, max: 5,
  palette: ['1a6b1a', '8db56e', '2e6fa3', 'f5d57c', 'd4b483']
  //          DenseF    DegradedF  Water    Agri/Sett  Bare
};

epochs.forEach(function(ep) {
  var composite  = buildComposite(ep);
  var classified = composite.classify(rf);

  // --- Export classified GeoTIFF ---
  Export.image.toDrive({
    image:       classified,
    description: 'RF_Classified_' + ep.year,
    folder:      'Rangamati_Deforestation',
    fileNamePrefix: 'RF_Classified_' + ep.year,
    region:      roi,
    scale:       30,
    crs:         'EPSG:32646',  // UTM Zone 46N — proper for Rangamati
    maxPixels:   1e10,
    fileFormat:  'GeoTIFF'
  });

  // --- Export forest binary mask (1=Forest, 0=Non-forest) ---
  var forestMask = classified.lte(2).selfMask();  // classes 1+2 = forest
  Export.image.toDrive({
    image:       forestMask.rename('forest'),
    description: 'Forest_Mask_' + ep.year,
    folder:      'Rangamati_Deforestation',
    fileNamePrefix: 'Forest_Mask_' + ep.year,
    region:      roi,
    scale:       30,
    crs:         'EPSG:32646',
    maxPixels:   1e10,
    fileFormat:  'GeoTIFF'
  });

  print('Export submitted: RF_Classified_' + ep.year);
});

// ================================================================
// PART 2: EXPORT 2023 RGB COMPOSITE (for basemap in publication figures)
// ================================================================
var composite2023vis = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2023-01-01', '2023-04-30')
  .map(maskAndScale)
  .map(harmonizeL8)
  .median()
  .clip(roi);

Export.image.toDrive({
  image:          composite2023vis.select(['Red','Green','Blue']).multiply(3.5).clamp(0,1),
  description:    'Composite_2023_TrueColor',
  folder:         'Rangamati_Deforestation',
  fileNamePrefix: 'Composite_2023_TrueColor',
  region:         roi,
  scale:          30,
  crs:            'EPSG:32646',
  maxPixels:      1e10,
  fileFormat:     'GeoTIFF'
});

// ================================================================
// PART 3: EXPORT NDVI CHANGE (1993 vs 2023) for visualization
// ================================================================
var comp1993 = buildComposite({
  year: 1993, sensor: 'L5',
  col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5
});

var ndviChange = composite2023.select('NDVI')
  .subtract(comp1993.select('NDVI'))
  .rename('NDVI_change_1993_2023')
  .clip(roi);

Export.image.toDrive({
  image:          ndviChange,
  description:    'NDVI_Change_1993_2023',
  folder:         'Rangamati_Deforestation',
  fileNamePrefix: 'NDVI_Change_1993_2023',
  region:         roi,
  scale:          30,
  crs:            'EPSG:32646',
  maxPixels:      1e10,
  fileFormat:     'GeoTIFF'
});

// ================================================================
// PART 4: ADD VISUAL LAYERS TO MAP (for checking before export)
// ================================================================
var composite2023_check = buildComposite({
  year: 2023, sensor: 'L8',
  col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8
});
var classified2023 = composite2023_check.classify(rf);
var classified1993 = buildComposite({
  year: 1993, sensor: 'L5',
  col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5
}).classify(rf);

Map.centerObject(roi, 10);
Map.addLayer(composite2023vis.select(['Red','Green','Blue']),
  {min: 0, max: 0.25}, '2023 True Color');
Map.addLayer(classified2023, classVis, '2023 Classified (RF)');
Map.addLayer(classified1993, classVis, '1993 Classified (RF)');
Map.addLayer(ndviChange, {min: -0.4, max: 0.4,
  palette: ['d73027','fc8d59','fee08b','d9ef8b','91cf60','1a9850']},
  'NDVI Change 1993→2023');

print('All exports submitted. Check Tasks tab (top-right) for progress.');
print('Files will appear in Google Drive > Rangamati_Deforestation folder.');
print('Total exports: ' + (6 + 6 + 1 + 1) + ' files (classified + forest masks + RGB + NDVI change)');
