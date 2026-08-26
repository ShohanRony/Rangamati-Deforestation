// ================================================================
// RANGAMATI DEFORESTATION MONITORING
// Script 07: Class Area Calculation for Olofsson Accuracy
// NO EXPORTS — Console output only
// ================================================================

// ---- STUDY AREA ----
var studyArea = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
).filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
var roi = studyArea.geometry();

// ---- PREPROCESSING ----
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

// ---- BUILD 2023 COMPOSITE ----
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

// ---- TRAINING DATA ----
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

samples = samples.map(function(f) {
  var c  = f.geometry().coordinates();
  var bx = ee.Number(c.get(0)).multiply(10).floor();
  var by = ee.Number(c.get(1)).multiply(10).floor();
  var block = bx.multiply(10000).add(by);
  return f.set('block_id', block).set('block_fold', block.abs().mod(10));
});

var trainSet = samples.filter(ee.Filter.lt('block_fold', 7));

// ---- TRAIN RF ----
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({
    features:        trainSet,
    classProperty:   'landcover',
    inputProperties: featureBands
  });

// ---- CLASSIFY 2023 ----
var classified2023 = composite2023.classify(rf);

// ---- COMPUTE AREAS ----
var areaImage = ee.Image.pixelArea()
  .divide(10000)
  .addBands(classified2023.rename('class'));

var classAreas = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName:  'class'
  }),
  geometry:  roi,
  scale:     30,
  maxPixels: 1e13
});

// ---- PRINT TO CONSOLE ----
var classNames = {
  1: 'Dense Forest',
  2: 'Degraded/Jhum',
  3: 'Water',
  4: 'Agriculture/Settlement',
  5: 'Bare Land'
};

ee.List(classAreas.get('groups')).evaluate(function(list) {
  print('=== 2023 RF Classification — Class Areas ===');
  var total = 0;
  list.forEach(function(item) {
    total += item['sum'];
    print('Class ' + item['class'] + 
          ' (' + classNames[item['class']] + '): ' +
          item['sum'].toFixed(2) + ' ha');
  });
  print('─────────────────────────────────────────');
  print('Total Study Area: ' + total.toFixed(2) + ' ha');
});

// Map preview (optional — remove if slow)
Map.centerObject(roi, 10);
Map.addLayer(classified2023, {
  min: 1, max: 5,
  palette: ['1a6b1a','8db56e','2e6fa3','f5d57c','d4b483']
}, '2023 Classified');
