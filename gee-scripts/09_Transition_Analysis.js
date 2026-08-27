// ================================================================
// RANGAMATI DEFORESTATION MONITORING
// Script 09: Land-Cover Change and Transition Analysis
//
// Computes pairwise transition matrices between consecutive epochs
// (1993-1998, 1998-2003, 2003-2008, 2008-2018, 2018-2023) plus the
// cumulative 1993-2023 comparison, using the same spatial-block
// classifier pipeline as Scripts 05/07/08. 2013 is excluded — it is
// legacy/exploratory and not part of the final six-epoch analysis.
//
// Area per transition is computed with an ee.Image.pixelArea()-weighted
// grouped sum, giving the true geodesic area of every pixel regardless
// of its position within the study area's projection (consistent with
// the method used in Script 07).
//
// Because intervals are unequal (5, 5, 5, 10, 5 years), raw hectare
// changes between periods are NOT directly comparable — this script
// reports both raw area change AND annualized rate (ha/year) for every
// transition category.
//
// Classes: 1=Dense Forest 2=Degraded/Jhum 3=Water 4=Agri/Settlement 5=Bare Land
// ================================================================

var studyArea = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
).filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
var roi = studyArea.geometry();

// ---- PREPROCESSING (identical to Scripts 05/07/08) ----
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

var srtm      = ee.Image('USGS/SRTMGL1_003').clip(roi);
var elevation = srtm.select('elevation');
var slope     = ee.Terrain.slope(srtm);
var featureBands = ['Blue','Green','Red','NIR','SWIR1','SWIR2',
                    'NDVI','NDWI','EVI','NBR','elevation','slope'];

var epochs = [
  {year: 1993, col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5},
  {year: 1998, col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5},
  {year: 2003, col: 'LANDSAT/LE07/C02/T1_L2', harmonize: harmonizeL7},
  {year: 2008, col: 'LANDSAT/LT05/C02/T1_L2', harmonize: harmonizeL5},
  {year: 2018, col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8},
  {year: 2023, col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8},
];

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

// ---- REBUILD CLASSIFIER (identical to Scripts 05/07/08) ----
var worldcover = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(roi);
var composite2023 = buildComposite({year: 2023, col: 'LANDSAT/LC08/C02/T1_L2', harmonize: harmonizeL8});
var ndvi2023 = composite2023.select('NDVI');

var treeMask = worldcover.eq(10);
var dense    = treeMask.and(ndvi2023.gte(0.55));
var degraded = treeMask.and(ndvi2023.gte(0.25)).and(ndvi2023.lt(0.55))
               .or(worldcover.eq(20)).or(worldcover.eq(30));

var classMap = ee.Image(0)
  .where(dense,                                                1)
  .where(degraded,                                             2)
  .where(worldcover.eq(80),                                    3)
  .where(worldcover.eq(40).or(worldcover.eq(50)),              4)
  .where(worldcover.eq(60),                                    5)
  .rename('landcover')
  .updateMask(worldcover.gt(0));
classMap = classMap.updateMask(classMap.neq(0));

var samples = classMap.addBands(composite2023).stratifiedSample({
  numPoints: 150, classBand: 'landcover', region: roi,
  scale: 30, seed: 42, geometries: true
});

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

var rf = ee.Classifier.smileRandomForest({numberOfTrees: 500, seed: 42})
  .train({features: trainSet, classProperty: 'landcover', inputProperties: featureBands});

// ---- CLASSIFY ALL SIX EPOCHS ----
var classified = {};
epochs.forEach(function(ep) {
  var composite = (ep.year === 2023) ? composite2023 : buildComposite(ep);
  classified[ep.year] = composite.classify(rf);
});

// ----------------------------------------------------------------
// TRANSITION MATRIX FUNCTION
// Area per transition class is computed via an ee.Image.pixelArea()-
// weighted grouped sum, robust to grid/projection mismatches.
// ----------------------------------------------------------------
var classNames = {1:'Dense Forest', 2:'Degraded/Jhum', 3:'Water', 4:'Agri/Settlement', 5:'Bare Land'};

function runTransition(yearFrom, yearTo, years, label) {
  var imgFrom = classified[yearFrom];
  var imgTo   = classified[yearTo];
  var code = imgFrom.multiply(10).add(imgTo).rename('code');

  var areaImg = ee.Image.pixelArea().divide(10000).rename('area_ha').addBands(code);

  var grouped = areaImg.reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'code'}),
    geometry: roi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 4
  });

  grouped.get('groups').evaluate(function(groups) {
    print('════════════════════════════════════════════════════════════');
    print(label + '  (' + yearFrom + ' -> ' + yearTo + ', ' + years + ' years)');
    print('════════════════════════════════════════════════════════════');

    var matrix = {}, rowTotal = {}, colTotal = {};
    [1,2,3,4,5].forEach(function(i) {
      matrix[i] = {1:0,2:0,3:0,4:0,5:0};
      rowTotal[i] = 0; colTotal[i] = 0;
    });

    var grandTotal = 0;
    groups.forEach(function(g) {
      var code_num = g.code;
      var from = Math.floor(code_num / 10);
      var to   = code_num % 10;
      if (from < 1 || from > 5 || to < 1 || to > 5) return; // skip nodata edge codes
      var ha = g.sum;
      matrix[from][to] += ha;
      rowTotal[from] += ha;
      colTotal[to]   += ha;
      grandTotal     += ha;
    });

    function pad(s, n) {
      s = String(s);
      while (s.length < n) { s = s + ' '; }
      return s;
    }
    var header = '  from\\to : ';
    [1,2,3,4,5].forEach(function(j) { header += pad(classNames[j], 16); });
    print(header);
    [1,2,3,4,5].forEach(function(i) {
      var row = '  ' + pad(classNames[i], 10) + ': ';
      [1,2,3,4,5].forEach(function(j) {
        row += pad(Math.round(matrix[i][j]), 16);
      });
      print(row);
    });

    print('--- Per-class summary (ha) ---');
    [1,2,3,4,5].forEach(function(i) {
      var persistence = matrix[i][i];
      var loss = rowTotal[i] - persistence;
      var gain = colTotal[i] - persistence;
      var net  = colTotal[i] - rowTotal[i];
      print('  ' + classNames[i] + ':  ' +
            yearFrom + ' area=' + Math.round(rowTotal[i]) + 'ha, ' +
            yearTo + ' area=' + Math.round(colTotal[i]) + 'ha, ' +
            'persistence=' + Math.round(persistence) + 'ha, ' +
            'loss=' + Math.round(loss) + 'ha, gain=' + Math.round(gain) + 'ha, ' +
            'net=' + (net >= 0 ? '+' : '') + Math.round(net) + 'ha ' +
            '(' + (net/years >= 0 ? '+' : '') + Math.round(net/years) + ' ha/yr)');
    });

    var forestToNonforest = 0, nonforestToForest = 0;
    [1,2].forEach(function(f) {
      [3,4,5].forEach(function(nf) { forestToNonforest += matrix[f][nf]; });
    });
    [3,4,5].forEach(function(nf) {
      [1,2].forEach(function(f) { nonforestToForest += matrix[nf][f]; });
    });
    print('--- Forest / Non-forest conversion ---');
    print('  Forest -> Non-forest (deforestation): ' + Math.round(forestToNonforest) + ' ha (' +
          (forestToNonforest/years).toFixed(1) + ' ha/yr)');
    print('  Non-forest -> Forest (afforestation/regrowth): ' + Math.round(nonforestToForest) + ' ha (' +
          (nonforestToForest/years).toFixed(1) + ' ha/yr)');
    print('  Net forest change: ' + Math.round(nonforestToForest - forestToNonforest) + ' ha (' +
          ((nonforestToForest - forestToNonforest)/years).toFixed(1) + ' ha/yr)');

    print('--- Dense Forest <-> Degraded/Jhum (within-forest change) ---');
    print('  Dense -> Degraded (degradation): ' + Math.round(matrix[1][2]) + ' ha (' +
          (matrix[1][2]/years).toFixed(1) + ' ha/yr)');
    print('  Degraded -> Dense (regrowth/reclassification): ' + Math.round(matrix[2][1]) + ' ha (' +
          (matrix[2][1]/years).toFixed(1) + ' ha/yr)');

    print('--- Water transitions (CAUTION: may reflect Kaptai reservoir level fluctuation, not deforestation) ---');
    print('  Water -> other: ' + Math.round(rowTotal[3] - matrix[3][3]) + ' ha');
    print('  Other -> Water: ' + Math.round(colTotal[3] - matrix[3][3]) + ' ha');

    print('  Total area accounted: ' + Math.round(grandTotal) + ' ha (expect ~576,534 ha)');
    print('');
  });
}

// ---- RUN ALL FIVE CONSECUTIVE-PERIOD TRANSITIONS ----
runTransition(1993, 1998, 5,  'Period 1');
runTransition(1998, 2003, 5,  'Period 2');
runTransition(2003, 2008, 5,  'Period 3');
runTransition(2008, 2018, 10, 'Period 4 (10-year interval)');
runTransition(2018, 2023, 5,  'Period 5');

// ---- CUMULATIVE 1993 -> 2023 (30-year net change, for headline figures) ----
runTransition(1993, 2023, 30, 'CUMULATIVE (full 30-year period)');
