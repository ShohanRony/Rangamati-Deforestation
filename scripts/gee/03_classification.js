// ============================================================================
// Rangamati Land-Cover Change, 1993-2023
// Script 03 - Training-Sample Collection Prototype
// Researcher: Shohinur Pervez Shohan, RMSTU
// ----------------------------------------------------------------------------
// STATUS: preliminary / exploratory. Not part of the production pipeline.
//   Retained for provenance. Script 04 does not import this script's
//   output; it performs its own independent stratifiedSample() call
//   (same seed, 42) against ee.Image('ESA/WorldCover/v200/2021') — an
//   explicitly versioned single-image asset — rather than this script's
//   ee.ImageCollection('ESA/WorldCover/v200') mosaic.
//
// Purpose:  Prototype stratified random sampling of training points from
//           ESA WorldCover 2021, seed=42.
// Inputs:   projects/crypto-hallway-405211/assets/BGD_adm2 (private asset,
//           see docs/REPRODUCIBILITY.md); ESA/WorldCover/v200 (public).
// Outputs:  In-session FeatureCollection + console prints. No Drive export.
// Depends:  none.
// ============================================================================


// ============================================================
// SECTION 1: STUDY AREA
// ============================================================

var studyArea = ee.FeatureCollection('projects/crypto-hallway-405211/assets/BGD_adm2')
  .filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));

Map.centerObject(studyArea, 9);
Map.addLayer(studyArea, {color: 'white'}, 'Rangamati Boundary');


// ============================================================
// SECTION 2: LANDSAT SCALING FUNCTIONS (Collection 2, Level 2)
// ============================================================

function scaleL5(img) {
  var optical = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true)
            .copyProperties(img, img.propertyNames());
}

function scaleL7(img) {
  var optical = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true)
            .copyProperties(img, img.propertyNames());
}

function scaleL8(img) {
  var optical = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true)
            .copyProperties(img, img.propertyNames());
}

// Roy et al. (2016): Landsat 8 to Landsat 5 spectral harmonization
function harmonizeL8toL5(img) {
  var slopes     = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]);
  var intercepts = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]);
  var harmonized = img
    .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'])
    .multiply(slopes).add(intercepts)
    .rename(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']);
  return img.addBands(harmonized, null, true)
            .copyProperties(img, img.propertyNames());
}


// ============================================================
// SECTION 3: SPECTRAL INDICES
// ============================================================

// For Landsat 5 and 7: B4=NIR, B3=Red, B2=Green, B1=Blue
function addIndicesL5(img) {
  var ndvi = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI');
  var ndwi = img.normalizedDifference(['SR_B2','SR_B4']).rename('NDWI');
  var evi  = img.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR':  img.select('SR_B4'),
      'RED':  img.select('SR_B3'),
      'BLUE': img.select('SR_B1')
    }).rename('EVI');
  return img.addBands([ndvi, ndwi, evi]);
}

// For Landsat 8 after harmonization: bands are renamed to L5 names
function addIndicesL8(img) {
  var ndvi = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI');
  var ndwi = img.normalizedDifference(['SR_B2','SR_B4']).rename('NDWI');
  var evi  = img.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR':  img.select('SR_B4'),
      'RED':  img.select('SR_B3'),
      'BLUE': img.select('SR_B1')
    }).rename('EVI');
  return img.addBands([ndvi, ndwi, evi]);
}


// ============================================================
// SECTION 4: COMPOSITE BUILDERS
// ============================================================

// Landsat 5: 1993, 1998, 2008 (Jan–Apr dry season)
function makeL5Composite(year) {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 30))
    .map(scaleL5)
    .map(addIndicesL5)
    .median()
    .clip(studyArea)
    .set('year', year);
}

// Landsat 7: 2003 only (before SLC failure in May 2003)
function makeL7Composite(year) {
  return ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year + '-01-01', year + '-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL7)
    .map(addIndicesL5)
    .median()
    .clip(studyArea)
    .set('year', year);
}

// Landsat 8: 2013, 2018, 2023
// Note: 2013 uses Apr–Jul because L8 launched Feb 2013
function makeL8Composite(year) {
  var start = (year === 2013) ? year + '-04-01' : year + '-01-01';
  var end   = (year === 2013) ? year + '-07-31' : year + '-04-30';
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUD_COVER', 30))
    .map(scaleL8)
    .map(harmonizeL8toL5)
    .map(addIndicesL8)
    .median()
    .clip(studyArea)
    .set('year', year);
}


// ============================================================
// SECTION 5: BUILD ALL 7 COMPOSITES
// ============================================================

var composite1993 = makeL5Composite(1993);
var composite1998 = makeL5Composite(1998);
var composite2003 = makeL7Composite(2003);
var composite2008 = makeL5Composite(2008);
var composite2013 = makeL8Composite(2013);
var composite2018 = makeL8Composite(2018);
var composite2023 = makeL8Composite(2023);


// ============================================================
// SECTION 6: DISPLAY 2023 AS VISUAL REFERENCE
// ============================================================

Map.addLayer(composite2023,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 0.0, max: 0.4},
  'False Color 2023 (NIR)');

Map.addLayer(composite2023,
  {bands: ['SR_B3','SR_B2','SR_B1'], min: 0.02, max: 0.25},
  'True Color 2023', false);

Map.addLayer(composite2023.select('NDVI'),
  {min: -0.1, max: 0.8, palette: ['red','yellow','green']},
  'NDVI 2023', false);


// ============================================================
// SECTION 7: TRAINING SAMPLE COLLECTION
// Source: ESA WorldCover 2021 (10m) + NDVI-based forest refinement
// ============================================================

// Load ESA WorldCover 2021
var worldcover = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .clip(studyArea);

// Use 2023 NDVI to separate dense vs degraded forest
// ESA classifies both as "Tree cover" (class 10) without density distinction
var ndvi2023 = composite2023.select('NDVI');
var esaTreeMask = worldcover.eq(10);

// Dense Forest: tree pixels with high NDVI (>= 0.55)
var denseMask    = esaTreeMask.and(ndvi2023.gte(0.55));

// Degraded Forest: tree pixels with moderate NDVI (0.25 – 0.55)
var degradedMask = esaTreeMask.and(ndvi2023.gt(0.25).and(ndvi2023.lt(0.55)));

// Build land cover class map
// Class 1 = Dense Forest
// Class 2 = Degraded Forest
// Class 3 = Water (Kaptai Lake + rivers)
// Class 4 = Agriculture / Settlement
// Class 5 = Jhum Cultivation (shrub + grassland)
// Class 6 = Bare Land / Other

var classMap = ee.Image(0)
  .where(denseMask,                                    1)
  .where(degradedMask,                                 2)
  .where(worldcover.eq(80),                            3)
  .where(worldcover.eq(40).or(worldcover.eq(50)),      4)
  .where(worldcover.eq(20).or(worldcover.eq(30)),      5)
  .where(worldcover.eq(60),                            6)
   .clip(studyArea)
  .rename('landcover');
  classMap = classMap.updateMask(classMap.neq(0)); // Exclude unclassified pixels

// Stratified random sampling: 100 points per class
// Fixed seed = 42 ensures reproducibility across runs
var trainingSamples = classMap.stratifiedSample({
  numPoints: 100,
  classBand: 'landcover',
  region:    studyArea,
  scale:     30,
  seed:      42,
  geometries: true
});

// Console output
print('Total training points:', trainingSamples.size());
print('Points per class:', trainingSamples.aggregate_histogram('landcover'));
print('Class legend: 1=Dense Forest | 2=Degraded Forest | 3=Water | 4=Agri/Settlement | 5=Jhum | 6=Bare Land');

// Display class map and training points
var classPalette = ['#1a5e1a','#8db36a','#1a3cff','#ffcc00','#ff8c00','#c8a882'];

Map.addLayer(classMap,
  {min: 1, max: 6, palette: classPalette},
  'Land Cover Class Map', false);

Map.addLayer(trainingSamples,
  {color: 'yellow'},
  'Training Points (stratified)');
