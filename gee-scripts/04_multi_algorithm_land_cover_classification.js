// ============================================================
// PROJECT: Deforestation Monitoring - Rangamati
// Script 04: Multi-Algorithm Land Cover Classification
// Researcher: Shohinur Pervez Shohan, RMSTU
// Classes: 5 (Dense Forest, Degraded Forest/Shrub, Water,
//             Agriculture/Settlement, Bare Land)
// Validation: Spatial-block holdout (0.1 deg blocks, ~30% of blocks held out)
//   NOTE (fixed): each 0.1-deg block is assigned WHOLLY to train or holdout
//   via a seeded random draw on the block, not on individual points, and not
//   via a coordinate-derived modulo (the previous block_id/mod(10) scheme
//   collapsed to a latitude-only split — see comments below).
// ============================================================

var studyArea = ee.FeatureCollection('projects/crypto-hallway-405211/assets/BGD_adm2')
  .filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
Map.centerObject(studyArea, 9);
Map.addLayer(studyArea, {color: 'white'}, 'Study area');

// ---------------------------------------------------------------------------
// Landsat preprocessing
// ---------------------------------------------------------------------------
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

// Roy et al. (2016): OLI -> ETM+ harmonization
function harmonizeL8(img) {
  var slopes     = ee.Image.constant([0.8850,0.9317,0.9372,0.8339,0.8639,0.9165]);
  var intercepts = ee.Image.constant([0.0183,0.0123,0.0123,0.0448,0.0306,0.0116]);
  var l8  = img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']);
  var out = l8.multiply(slopes).add(intercepts)
              .rename(['Blue','Green','Red','NIR','SWIR1','SWIR2']);
  return out.copyProperties(img, img.propertyNames());
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

function makeL57Composite(year) {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-05-01')
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(maskAndScale)
    .map(function(img) {
      return img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'])
                .rename(['Blue','Green','Red','NIR','SWIR1','SWIR2']);
    })
    .map(addIndices)
    .median().clip(studyArea).set('year', year);
}

function makeL8Composite(year) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-05-01')
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(maskAndScale).map(harmonizeL8).map(addIndices)
    .median().clip(studyArea).set('year', year);
}

var srtm    = ee.Image('USGS/SRTMGL1_003').clip(studyArea);
var terrain = srtm.rename('elevation')
                  .addBands(ee.Terrain.slope(srtm).rename('slope'));
function withTerrain(img) { return img.addBands(terrain); }

var composite1993 = withTerrain(makeL57Composite(1993));
var composite1998 = withTerrain(makeL57Composite(1998));
var composite2003 = withTerrain(makeL57Composite(2003));
var composite2008 = withTerrain(makeL57Composite(2008));
var composite2013 = withTerrain(makeL8Composite(2013));
var composite2018 = withTerrain(makeL8Composite(2018));
var composite2023 = withTerrain(makeL8Composite(2023));

var featureBands = ['Blue','Green','Red','NIR','SWIR1','SWIR2',
                    'NDVI','NDWI','EVI','NBR','elevation','slope'];

// ---------------------------------------------------------------------------
// Class map — 5 classes
// 1=Dense Forest  2=Degraded Forest/Shrub/Grass  3=Water
// 4=Agriculture/Settlement  5=Bare Land
// ---------------------------------------------------------------------------
var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(studyArea);
var ndvi2023   = composite2023.select('NDVI');
var tree       = worldcover.eq(10);

var dense    = tree.and(ndvi2023.gte(0.55));
var degraded = tree.and(ndvi2023.gte(0.25).and(ndvi2023.lt(0.55)))
                   .or(worldcover.eq(20)).or(worldcover.eq(30));

var rawMap = ee.Image(0)
  .where(dense,                                   1)
  .where(degraded,                                2)
  .where(worldcover.eq(80),                       3)
  .where(worldcover.eq(40).or(worldcover.eq(50)), 4)
  .where(worldcover.eq(60),                       5);

var classMap = rawMap.updateMask(rawMap.neq(0)).rename('landcover');

// ---------------------------------------------------------------------------
// Stratified random sampling (150 pts/class, seed=42)
// ---------------------------------------------------------------------------
var samples = classMap.stratifiedSample({
  numPoints: 150, classBand: 'landcover', region: studyArea,
  scale: 30, seed: 42, geometries: true, tileScale: 4
});

print('Total samples:', samples.size());
print('Per-class distribution:', samples.aggregate_histogram('landcover'));
print('Class legend: 1=Dense Forest | 2=Degraded/Shrub | 3=Water | 4=Agri/Settlement | 5=Bare Land');

// ---------------------------------------------------------------------------
// Spatial block assignment — FIXED
//
// Bug in the previous version: block_id was computed as
//   block = bx * 10000 + by
// and the fold was block_id.mod(10). Because 10000 is divisible by 10,
// (bx*10000 + by) mod 10 always equals (by mod 10) — the bx (longitude)
// term contributes nothing to the modulo. The "10 spatial blocks" were in
// fact 10 latitude stripes; longitude played no role in the split.
//
// Fix: build a 0.1-deg block id from (bx, by) as before (this id is now
// used ONLY as a block identifier, never fed into a modulo), then assign
// each DISTINCT block — not each point — to train or holdout via a single
// seeded random draw per block. This guarantees every point inside a given
// 0.1-deg cell goes to the same partition, and the partition itself is a
// genuine 2-D spatial split rather than a coordinate-derived pattern.
// ---------------------------------------------------------------------------
samples = samples.map(function(f) {
  var c   = f.geometry().coordinates();
  var bx  = ee.Number(c.get(0)).multiply(10).floor();
  var by  = ee.Number(c.get(1)).multiply(10).floor();
  var bid = bx.format('%d').cat('_').cat(by.format('%d'));
  return f.set('block_id', bid);
});

// One random draw per DISTINCT block (seed=42), ~30% of blocks -> holdout.
var uniqueBlocks = samples.distinct('block_id').randomColumn('rand', 42);
var blockFolds = uniqueBlocks.map(function(b) {
  var isHoldout = ee.Number(b.get('rand')).gte(0.7);
  return b.set('block_fold', ee.Algorithms.If(isHoldout, 1, 0));
});

print('Number of distinct 0.1-deg blocks:', uniqueBlocks.size());
print('Blocks assigned to holdout (fold=1):',
  blockFolds.filter(ee.Filter.eq('block_fold', 1)).size());
print('Blocks assigned to training (fold=0):',
  blockFolds.filter(ee.Filter.eq('block_fold', 0)).size());

// Join each point back to its block's fold assignment.
var joinFilter = ee.Filter.equals({leftField: 'block_id', rightField: 'block_id'});
var joined = ee.Join.saveFirst('blockMatch').apply(samples, blockFolds, joinFilter);
samples = joined.map(function(f) {
  var match = ee.Feature(f.get('blockMatch'));
  return f.set('block_fold', match.get('block_fold'));
});

var sampled = composite2023.select(featureBands).sampleRegions({
  collection: samples,
  properties: ['landcover','block_id','block_fold'],
  scale: 30, tileScale: 8, geometries: true
}).filter(ee.Filter.notNull(featureBands));

var trainSet = sampled.filter(ee.Filter.eq('block_fold', 0));
var testSet  = sampled.filter(ee.Filter.eq('block_fold', 1));

print('Training samples:', trainSet.size());
print('Validation (holdout) samples:', testSet.size());
print('Training class distribution:', trainSet.aggregate_histogram('landcover'));
print('Validation class distribution:', testSet.aggregate_histogram('landcover'));

// ---------------------------------------------------------------------------
// SVM z-score standardisation (fitted on training set only)
// ---------------------------------------------------------------------------
var zBands = featureBands.map(function(b) { return b + '_z'; });
var means  = featureBands.map(function(b) {
  return ee.Number(trainSet.aggregate_mean(b));
});
var sds = featureBands.map(function(b) {
  return ee.Number(trainSet.aggregate_total_sd(b)).max(0.000001);
});

function standardizeFC(fc) {
  return fc.map(function(f) {
    var vals = featureBands.map(function(b, i) {
      return ee.Number(f.get(b)).subtract(means[i]).divide(sds[i]);
    });
    return f.set(ee.Dictionary.fromLists(zBands, vals));
  });
}

function standardizeImg(img) {
  var bands = featureBands.map(function(b, i) {
    return img.select(b)
      .subtract(ee.Image.constant(means[i]))
      .divide(ee.Image.constant(sds[i]))
      .rename(b + '_z');
  });
  return ee.Image.cat(bands);
}

var trainZ = standardizeFC(trainSet);
var testZ  = standardizeFC(testSet);

// ---------------------------------------------------------------------------
// Train classifiers
// ---------------------------------------------------------------------------
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({features: trainSet, classProperty: 'landcover',
          inputProperties: featureBands});

var cart = ee.Classifier.smileCart(50, 3)
  .train({features: trainSet, classProperty: 'landcover',
          inputProperties: featureBands});

var svm = ee.Classifier.libsvm({kernelType: 'RBF', gamma: 0.05, cost: 10})
  .train({features: trainZ, classProperty: 'landcover',
          inputProperties: zBands});

// ---------------------------------------------------------------------------
// Spatial-block holdout validation
// ---------------------------------------------------------------------------
var rfMatrix   = testSet.classify(rf).errorMatrix('landcover', 'classification');
var cartMatrix = testSet.classify(cart).errorMatrix('landcover', 'classification');
var svmMatrix  = testZ.classify(svm).errorMatrix('landcover', 'classification');

print('=== SPATIAL-BLOCK HOLDOUT VALIDATION ===');
print('RF   OA / Kappa:', rfMatrix.accuracy(),   rfMatrix.kappa());
print('CART OA / Kappa:', cartMatrix.accuracy(), cartMatrix.kappa());
print('SVM  OA / Kappa:', svmMatrix.accuracy(),  svmMatrix.kappa());
print('RF confusion matrix:',   rfMatrix);
print('CART confusion matrix:', cartMatrix);
print('SVM confusion matrix:',  svmMatrix);

// ---------------------------------------------------------------------------
// Classify 2023 + visualize
// ---------------------------------------------------------------------------
var rfClassified2023 = composite2023.select(featureBands).classify(rf);

var palette  = ['#1a5e1a','#8db36a','#1a3cff','#ffcc00','#c8a882'];
var visClass = {min: 1, max: 5, palette: palette};

Map.addLayer(composite2023.select('NDVI'),
  {min:0, max:1, palette:['brown','yellow','green']}, 'NDVI 2023', false);
Map.addLayer(rfClassified2023, visClass, 'RF 2023');
Map.addLayer(composite2023.select(featureBands).classify(cart), visClass, 'CART 2023', false);
Map.addLayer(standardizeImg(composite2023).classify(svm), visClass, 'SVM 2023', false);

// ---------------------------------------------------------------------------
// Apply RF to full time series + forest area by year
// ---------------------------------------------------------------------------
var composites = [composite1993, composite1998, composite2003,
                  composite2008, composite2013, composite2018, composite2023];
var years      = [1993, 1998, 2003, 2008, 2013, 2018, 2023];

var classified = composites.map(function(img) {
  return img.select(featureBands).classify(rf).set('year', img.get('year'));
});

var areaFeatures = years.map(function(year, i) {
  var forest = classified[i].eq(1).or(classified[i].eq(2));
  var area   = forest.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: studyArea.geometry(),
    scale: 30, maxPixels: 1e10, tileScale: 4
  });
  return ee.Feature(null, {
    year: year,
    forest_area_km2: ee.Number(area.get('classification')).divide(1e6)
  });
});

print('Forest area by year:', ee.FeatureCollection(areaFeatures));

// ---------------------------------------------------------------------------
// EXPORT: Confusion matrices + class areas (for Olofsson analysis)
// Filenames suffixed "_v2" so they don't overwrite the earlier
// (buggy-split) exports — keep both, but use only the _v2 files in the paper.
// ---------------------------------------------------------------------------
Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(null, {
    'classifier': 'RF',
    'OA':    rfMatrix.accuracy(),
    'kappa': rfMatrix.kappa(),
    'matrix': rfMatrix.array()
  })]),
  description: 'RF_ConfusionMatrix_Spatial_Block_v2',
  fileFormat: 'CSV', folder: 'Rangamati_Deforestation'
});

Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(null, {
    'classifier': 'CART',
    'OA':    cartMatrix.accuracy(),
    'kappa': cartMatrix.kappa(),
    'matrix': cartMatrix.array()
  })]),
  description: 'CART_ConfusionMatrix_Spatial_Block_v2',
  fileFormat: 'CSV', folder: 'Rangamati_Deforestation'
});

Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(null, {
    'classifier': 'SVM',
    'OA':    svmMatrix.accuracy(),
    'kappa': svmMatrix.kappa(),
    'matrix': svmMatrix.array()
  })]),
  description: 'SVM_ConfusionMatrix_Spatial_Block_v2',
  fileFormat: 'CSV', folder: 'Rangamati_Deforestation'
});

var classAreas2023 = ee.List([1,2,3,4,5]).map(function(c) {
  var area = rfClassified2023.eq(ee.Number(c))
    .multiply(ee.Image.pixelArea())
    .reduceRegion({
      reducer: ee.Reducer.sum(), geometry: studyArea.geometry(),
      scale: 30, maxPixels: 1e10, tileScale: 4
    });
  return ee.Feature(null, {
    'class':    c,
    'area_km2': ee.Number(area.get('classification')).divide(1e6)
  });
});

Export.table.toDrive({
  collection: ee.FeatureCollection(classAreas2023),
  description: 'RF_MappedArea_PerClass_2023_v2',
  fileFormat: 'CSV', folder: 'Rangamati_Deforestation'
});
