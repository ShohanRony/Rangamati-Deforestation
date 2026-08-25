// ============================================================
// PROJECT: Deforestation Monitoring - Rangamati
// Script 04: Multi-Algorithm Land Cover Classification
// Researcher: Shohinur Pervez Shohan, RMSTU
// Classifiers: Random Forest (500 trees), CART, SVM (RBF kernel)
// Validation:  Spatial-block holdout (0.1 deg blocks, 30% test split)
// Note: WorldCover/NDVI labels are provisional reference labels,
//       not field-verified truth. Audit against historical imagery
//       before citing results in a publication.
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
  // QA_PIXEL bits: 1=dilated cloud, 3=cloud, 4=cloud shadow, 5=snow.
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

// Roy et al. (2016) coefficients: convert OLI reflectance to ETM+-like
// reflectance before using the same feature names across sensors.
function harmonizeL8(img) {
  var slopes = ee.Image.constant([0.8850, 0.9317, 0.9372, 0.8339, 0.8639, 0.9165]);
  var intercepts = ee.Image.constant([0.0183, 0.0123, 0.0123, 0.0448, 0.0306, 0.0116]);
  var l8 = img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']);
  var out = l8.multiply(slopes).add(intercepts)
    .rename(['Blue','Green','Red','NIR','SWIR1','SWIR2']);
  return out.copyProperties(img, img.propertyNames());
}

function addIndices(img) {
  var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var ndwi = img.normalizedDifference(['Green', 'NIR']).rename('NDWI');
  var evi = img.expression(
    '2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))', {
      nir: img.select('NIR'), red: img.select('Red'), blue: img.select('Blue')
    }).rename('EVI');
  var nbr = img.normalizedDifference(['NIR', 'SWIR2']).rename('NBR');
  return img.addBands([ndvi, ndwi, evi, nbr]);
}

function makeL57Composite(year) {
  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .filterBounds(studyArea).filterDate(year + '-01-01', year + '-05-01')
    .filter(ee.Filter.lt('CLOUD_COVER', 70)).map(maskAndScale)
    .map(function(img) {
      return img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'])
        .rename(['Blue','Green','Red','NIR','SWIR1','SWIR2']);
    }).map(addIndices);
  return col.median().clip(studyArea).set('year', year);
}

function makeL8Composite(year) {
  var col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-12-31')  // full year for 2013
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(maskAndScale).map(harmonizeL8).map(addIndices);
  return col.median().clip(studyArea).set('year', year);
}

var srtm = ee.Image('USGS/SRTMGL1_003').clip(studyArea);
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
// Provisional labels (must be audited against historical imagery)
// ---------------------------------------------------------------------------
var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(studyArea);
var ndvi2023 = composite2023.select('NDVI');
var tree = worldcover.eq(10);
var dense = tree.and(ndvi2023.gte(0.55));
var degraded = tree.and(ndvi2023.gte(0.25).and(ndvi2023.lt(0.55)))
  .or(worldcover.eq(20)).or(worldcover.eq(30));

// 1=dense forest, 2=degraded forest/shrub/grass, 3=water,
// 4=agriculture/built-up, 5=bare/sparse vegetation.
var classMap = ee.Image(0)
  .where(dense, 1).where(degraded, 2)
  .where(worldcover.eq(80), 3)
  .where(worldcover.eq(40).or(worldcover.eq(50)), 4)
  .where(worldcover.eq(60), 5)
  .updateMask(ee.Image(0).where(dense, 1).where(degraded, 2)
    .where(worldcover.eq(80), 3)
    .where(worldcover.eq(40).or(worldcover.eq(50)), 4)
    .where(worldcover.eq(60), 5).neq(0))
  .rename('landcover');

var samples = classMap.stratifiedSample({
  numPoints: 150, classBand: 'landcover', region: studyArea,
  scale: 30, seed: 42, geometries: true, tileScale: 4
});

// Assign a 0.1-degree spatial block. The block, not the pixel, is the split
// unit, so nearby pixels cannot appear in both train and validation sets.
samples = samples.map(function(f) {
  var c = f.geometry().coordinates();
  var bx = ee.Number(c.get(0)).multiply(10).floor();
  var by = ee.Number(c.get(1)).multiply(10).floor();
  var block = bx.multiply(10000).add(by);
  return f.set('block_id', block).set('block_fold', block.abs().mod(10));
});
var sampled = composite2023.select(featureBands).sampleRegions({
  collection: samples, properties: ['landcover','block_id','block_fold'],
  scale: 30, tileScale: 8, geometries: true
}).filter(ee.Filter.notNull(featureBands));
var trainSet = sampled.filter(ee.Filter.lt('block_fold', 7));
var testSet = sampled.filter(ee.Filter.gte('block_fold', 7));
print('Spatial-block training samples:', trainSet.size());
print('Spatial-block validation samples:', testSet.size());
print('Class distribution:', sampled.aggregate_histogram('landcover'));
print('Training class distribution:', trainSet.aggregate_histogram('landcover'));
print('Validation class distribution:', testSet.aggregate_histogram('landcover'));

// ---------------------------------------------------------------------------
// SVM standardisation (parameters are fitted on training blocks only)
// ---------------------------------------------------------------------------
var zBands = featureBands.map(function(b) { return b + '_z'; });
var means = featureBands.map(function(b) { return ee.Number(trainSet.aggregate_mean(b)); });
var sds = featureBands.map(function(b) {
  return ee.Number(trainSet.aggregate_total_sd(b)).max(0.000001);
});
function standardizeFeatureCollection(fc) {
  return fc.map(function(f) {
    var vals = featureBands.map(function(b, i) {
      return ee.Number(f.get(b)).subtract(means[i]).divide(sds[i]);
    });
    return f.set(ee.Dictionary.fromLists(zBands, vals));
  });
}
function standardizedImage(img) {
  var bands = featureBands.map(function(b, i) {
    return img.select(b).subtract(ee.Image.constant(means[i]))
      .divide(ee.Image.constant(sds[i])).rename(b + '_z');
  });
  return ee.Image.cat(bands);
}
var trainZ = standardizeFeatureCollection(trainSet);
var testZ = standardizeFeatureCollection(testSet);

// ---------------------------------------------------------------------------
// Classifiers
// ---------------------------------------------------------------------------
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({features: trainSet, classProperty: 'landcover', inputProperties: featureBands});
var cart = ee.Classifier.smileCart(50, 3)
  .train({features: trainSet, classProperty: 'landcover', inputProperties: featureBands});
var svm = ee.Classifier.libsvm({kernelType: 'RBF', gamma: 0.05, cost: 10})
  .train({features: trainZ, classProperty: 'landcover', inputProperties: zBands});

var rfMatrix = testSet.classify(rf).errorMatrix('landcover', 'classification');
var cartMatrix = testSet.classify(cart).errorMatrix('landcover', 'classification');
var svmMatrix = testZ.classify(svm).errorMatrix('landcover', 'classification');
print('SPATIAL-BLOCK VALIDATION (not random-pixel CV)');
print('RF accuracy / kappa:', rfMatrix.accuracy(), rfMatrix.kappa());
print('CART accuracy / kappa:', cartMatrix.accuracy(), cartMatrix.kappa());
print('SVM accuracy / kappa:', svmMatrix.accuracy(), svmMatrix.kappa());
print('RF confusion matrix:', rfMatrix);
print('CART confusion matrix:', cartMatrix);
print('SVM confusion matrix:', svmMatrix);

var palette = ['#1a5e1a','#8db36a','#1a3cff','#ffcc00','#c8a882'];
var visClass = {min: 1, max: 5, palette: palette};
Map.addLayer(composite2023.select('NDVI'), {min: 0, max: 1, palette: ['brown','yellow','green']}, '2023 NDVI', false);
Map.addLayer(composite2023.select(featureBands).classify(rf), visClass, 'RF 2023');
Map.addLayer(composite2023.select(featureBands).classify(cart), visClass, 'CART 2023', false);
Map.addLayer(standardizedImage(composite2023).classify(svm), visClass, 'SVM 2023', false);

// Apply RF to the time series. Historical reference-label auditing is still
// required before interpreting these as validated historical maps.
var composites = [composite1993, composite1998, composite2003, composite2008,
                  composite2013, composite2018, composite2023];
var years = [1993, 1998, 2003, 2008, 2013, 2018, 2023];
var classified = composites.map(function(img) {
  return img.select(featureBands).classify(rf).set('year', img.get('year'));
});
var areaFeatures = years.map(function(year, i) {
  var forest = classified[i].eq(1).or(classified[i].eq(2));
  var area = forest.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: studyArea, scale: 30,
    maxPixels: 1e10, tileScale: 4
  });
  return ee.Feature(null, {year: year,
    forest_area_km2: ee.Number(area.get('classification')).divide(1e6)});
});
print('Forest area by year (provisional):', ee.FeatureCollection(areaFeatures))
