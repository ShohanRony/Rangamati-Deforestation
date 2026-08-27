// ============================================================
// SCRIPT 06a — Independent Visual Validation Points
// Complete script — no dependency on Script 04 session
// ============================================================

var studyArea = ee.FeatureCollection('projects/crypto-hallway-405211/assets/BGD_adm2')
  .filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
Map.centerObject(studyArea, 9);

// ── Preprocessing ─────────────────────────────────────────────
function maskAndScale(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1<<1).eq(0)
    .and(qa.bitwiseAnd(1<<3).eq(0))
    .and(qa.bitwiseAnd(1<<4).eq(0))
    .and(qa.bitwiseAnd(1<<5).eq(0));
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
    '2.5*((nir-red)/(nir+6*red-7.5*blue+1))',
    {nir:img.select('NIR'),red:img.select('Red'),blue:img.select('Blue')}
  ).rename('EVI');
  var nbr  = img.normalizedDifference(['NIR','SWIR2']).rename('NBR');
  return img.addBands([ndvi, ndwi, evi, nbr]);
}

var srtm    = ee.Image('USGS/SRTMGL1_003').clip(studyArea);
var terrain = srtm.rename('elevation')
  .addBands(ee.Terrain.slope(srtm).rename('slope'));

function makeL8Composite(year) {
  var col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year+'-01-01', year+'-05-01')
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(maskAndScale).map(harmonizeL8).map(addIndices);
  return col.median().clip(studyArea)
    .addBands(terrain).set('year', year);
}

var featureBands = ['Blue','Green','Red','NIR','SWIR1','SWIR2',
                    'NDVI','NDWI','EVI','NBR','elevation','slope'];

var composite2023 = makeL8Composite(2023);

// ── Training labels (same as Script 04) ───────────────────────
var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(studyArea);
var ndvi2023   = composite2023.select('NDVI');
var tree       = worldcover.eq(10);
var dense      = tree.and(ndvi2023.gte(0.55));
var degraded   = tree.and(ndvi2023.gte(0.25).and(ndvi2023.lt(0.55)))
                     .or(worldcover.eq(20)).or(worldcover.eq(30));

var rawMap = ee.Image(0)
  .where(dense, 1).where(degraded, 2)
  .where(worldcover.eq(80), 3)
  .where(worldcover.eq(40).or(worldcover.eq(50)), 4)
  .where(worldcover.eq(60), 5);
var classMap = rawMap.updateMask(rawMap.neq(0)).rename('landcover');

// ── Spatial-block training samples (same as Script 04) ────────
var samples = classMap.stratifiedSample({
  numPoints: 150, classBand: 'landcover', region: studyArea,
  scale: 30, seed: 42, geometries: true, tileScale: 4
});
samples = samples.map(function(f) {
  var c     = f.geometry().coordinates();
  var bx    = ee.Number(c.get(0)).multiply(10).floor();
  var by    = ee.Number(c.get(1)).multiply(10).floor();
  var block = bx.multiply(10000).add(by);
  return f.set('block_fold', block.abs().mod(10));
});
var sampled = composite2023.select(featureBands).sampleRegions({
  collection: samples, properties: ['landcover','block_fold'],
  scale: 30, tileScale: 8, geometries: true
}).filter(ee.Filter.notNull(featureBands));
var trainSet = sampled.filter(ee.Filter.lt('block_fold', 7));

// ── Train RF (same as Script 04) ──────────────────────────────
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({features: trainSet, classProperty: 'landcover',
          inputProperties: featureBands});

// ── Classify 2023 ─────────────────────────────────────────────
var classified2023 = composite2023.select(featureBands).classify(rf);

// ── Generate 200 independent validation points ─────────────────
// seed=99 — completely different from Script 04 seed=42
var validationPoints = classified2023.stratifiedSample({
  numPoints : 40,
  classBand : 'classification',
  region    : studyArea,
  scale     : 30,
  seed      : 99,
  geometries: true,
  tileScale : 4
});

var classNames = ee.Dictionary({
  '1': 'Dense Forest',
  '2': 'Degraded Forest / Jhum',
  '3': 'Water',
  '4': 'Agriculture / Settlement',
  '5': 'Bare Land'
});

validationPoints = validationPoints.map(function(f) {
  var cls = ee.Number(f.get('classification')).toInt().format('%d');
  return f
    .set('point_id',      f.id())
    .set('rf_class',      f.get('classification'))
    .set('rf_class_name', classNames.get(cls))
    .set('manual_class',  '')
    .set('manual_name',   '')
    .set('confidence',    '')
    .set('notes',         '');
});

print('Total points:', validationPoints.size());
print('Class distribution:', validationPoints.aggregate_histogram('rf_class'));

// ── Visualise ─────────────────────────────────────────────────
var palette = ['#1a5e1a','#8db36a','#1a3cff','#ffcc00','#c8a882'];
Map.addLayer(classified2023, {min:1, max:5, palette:palette}, 'RF 2023');
Map.addLayer(validationPoints, {color: 'red'}, 'Validation Points (200)');

// ── Export CSV + KML ──────────────────────────────────────────
Export.table.toDrive({
  collection : validationPoints,
  description: 'ValidationPoints_200_Script06a',
  folder     : 'Rangamati_Deforestation',
  fileFormat : 'CSV'
});

Export.table.toDrive({
  collection : validationPoints,
  description: 'ValidationPoints_200_KML',
  folder     : 'Rangamati_Deforestation',
  fileFormat : 'KML'
});

print('✓ Run both exports from Tasks tab.');
print('Then open KML in Google Earth Pro and fill manual_class in CSV.');
