// ============================================================================
// Rangamati Land-Cover Change, 1993-2023
// Script 02 - Early Multi-Epoch Compositing Pipeline (draft)
// Researcher: Shohinur Pervez Shohan, RMSTU
// ----------------------------------------------------------------------------
// STATUS: preliminary / exploratory. Not part of the production pipeline.
//   Retained for provenance. Differs from the production pipeline
//   (Script 04 onward) in three respects, noted here so results are not
//   mixed up across script versions:
//     1. Uses the Roy et al. (2016) Table 2 OLS coefficients in the
//        ETM+ -> OLI direction; the production scripts (04, 05, 07, 08, 09)
//        use the OLI -> ETM+ direction with different coefficient values.
//     2. Includes a 2013 Landsat 8 epoch; the production analysis uses six
//        epochs only (1993, 1998, 2003, 2008, 2018, 2023) and excludes 2013.
//     3. Produces a raw pixel-wise NDVI-difference "change" layer; the
//        production change analysis (Script 09) instead differences
//        classified maps and reports a full transition matrix.
//
// Purpose:  Build multi-year Landsat composites with spectral indices and
//           preview a naive NDVI-difference change layer.
// Inputs:   projects/crypto-hallway-405211/assets/BGD_adm2 (private asset,
//           see docs/REPRODUCIBILITY.md); Landsat 5/7/8 C02 T1_L2 collections.
// Outputs:  Map layers + console prints only. No exports.
// Depends:  none.
// ============================================================================

var bgd = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
);
var studyArea = bgd
  .filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'))
  .geometry();

// ================================================
// FUNCTIONS
// ================================================

function scaleL5(img) {
  var optical = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true).copyProperties(img, img.propertyNames());
}
function scaleL7(img) {
  var optical = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  return img.addBands(optical, null, true).copyProperties(img, img.propertyNames());
}
function scaleL8(img) {
  var optical = img.select('SR_B.').multiply(0.0000275).add(-0.2);
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

// L7 bands: B1=Blue, B2=Green, B3=Red, B4=NIR, B5=SWIR1, B7=SWIR2 (same as L5)
function addIndicesL5(img) {
  var ndvi = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI');
  var ndwi = img.normalizedDifference(['SR_B2','SR_B4']).rename('NDWI');
  var evi  = img.expression('2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
    {'NIR':img.select('SR_B4'),'RED':img.select('SR_B3'),'BLUE':img.select('SR_B1')}).rename('EVI');
  return img.addBands([ndvi, ndwi, evi]);
}
function addIndicesL8(img) {
  // After harmonization, L8 bands renamed to L5 names
  var ndvi = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI');
  var ndwi = img.normalizedDifference(['SR_B2','SR_B4']).rename('NDWI');
  var evi  = img.expression('2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
    {'NIR':img.select('SR_B4'),'RED':img.select('SR_B3'),'BLUE':img.select('SR_B1')}).rename('EVI');
  return img.addBands([ndvi, ndwi, evi]);
}

// ================================================
// BUILD COMPOSITES
// ================================================

// L5: 1993, 1998, 2008
function makeL5Composite(year) {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year+'-01-01', year+'-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL5).map(addIndicesL5)
    .median().clip(studyArea).set('year', year);
}

// L7: 2003 (before SLC failure in May 2003 — safe to use Jan-Apr)
function makeL7Composite(year) {
  return ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(year+'-01-01', year+'-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL7).map(addIndicesL5)  // L7 same band structure as L5
    .median().clip(studyArea).set('year', year);
}

// L8: 2013 (launched Feb 2013, data from Apr), 2018, 2023
function makeL8Composite(year) {
  var start = (year === 2013) ? year+'-04-01' : year+'-01-01';
  var end   = (year === 2013) ? year+'-07-31' : year+'-04-30';
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUD_COVER', 50))
    .map(scaleL8).map(harmonizeL8toL5).map(addIndicesL8)
    .median().clip(studyArea).set('year', year);
}

var comp1993 = makeL5Composite(1993);
var comp1998 = makeL5Composite(1998);
var comp2003 = makeL7Composite(2003);  // L7
var comp2008 = makeL5Composite(2008);
var comp2013 = makeL8Composite(2013);  // L8 (Apr-Jul)
var comp2018 = makeL8Composite(2018);
var comp2023 = makeL8Composite(2023);

// ================================================
// VERIFY
// ================================================
print('=== IMAGE COUNT ===');
print('L5-1993:', ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(studyArea).filterDate('1993-01-01','1993-04-30').filter(ee.Filter.lt('CLOUD_COVER',50)).size());
print('L5-1998:', ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(studyArea).filterDate('1998-01-01','1998-04-30').filter(ee.Filter.lt('CLOUD_COVER',50)).size());
print('L7-2003:', ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(studyArea).filterDate('2003-01-01','2003-04-30').filter(ee.Filter.lt('CLOUD_COVER',50)).size());
print('L5-2008:', ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(studyArea).filterDate('2008-01-01','2008-04-30').filter(ee.Filter.lt('CLOUD_COVER',50)).size());
print('L8-2013:', ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(studyArea).filterDate('2013-04-01','2013-07-31').filter(ee.Filter.lt('CLOUD_COVER',50)).size());
print('L8-2018:', ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(studyArea).filterDate('2018-01-01','2018-04-30').filter(ee.Filter.lt('CLOUD_COVER',50)).size());
print('L8-2023:', ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(studyArea).filterDate('2023-01-01','2023-04-30').filter(ee.Filter.lt('CLOUD_COVER',50)).size());

// ================================================
// VISUALIZATION
// ================================================
var ndviViz = {min:-0.1, max:0.8, palette:['brown','yellow','lightgreen','green','darkgreen']};

Map.centerObject(studyArea, 10);
Map.addLayer(bgd.filter(ee.Filter.eq('NAME_2','Parbattya Chattagram')), {color:'red'}, 'Rangamati Boundary');

Map.addLayer(comp1993.select('NDVI'), ndviViz, 'NDVI 1993', false);
Map.addLayer(comp1998.select('NDVI'), ndviViz, 'NDVI 1998', false);
Map.addLayer(comp2003.select('NDVI'), ndviViz, 'NDVI 2003', false);
Map.addLayer(comp2008.select('NDVI'), ndviViz, 'NDVI 2008', false);
Map.addLayer(comp2013.select('NDVI'), ndviViz, 'NDVI 2013', false);
Map.addLayer(comp2018.select('NDVI'), ndviViz, 'NDVI 2018', false);
Map.addLayer(comp2023.select('NDVI'), ndviViz, 'NDVI 2023', true);

Map.addLayer(comp1993,
  {bands:['SR_B3','SR_B2','SR_B1'], min:0.02, max:0.25},
  'True Color 1993', false);
Map.addLayer(comp2023,
  {bands:['SR_B3','SR_B2','SR_B1'], min:0.02, max:0.25},
  'True Color 2023', false);

print('=== PIPELINE COMPLETE ===');
print('Sensors: L5(1993,1998,2008) | L7(2003) | L8(2013,2018,2023)');

// === NDVI Change Detection (1993 vs 2023) ===
var ndviChange = comp2023.select('NDVI')
  .subtract(comp1993.select('NDVI'))
  .rename('NDVI_Change');

Map.addLayer(ndviChange, {
  min: -0.5, max: 0.5,
  palette: ['red','orange','white','lightgreen','darkgreen']
}, 'NDVI Change 1993→2023', true);

print('NDVI Change layer added');
print('Red = forest loss, Green = forest gain, White = no change');

// NDVI statistics check
print('1993 NDVI mean:', comp1993.select('NDVI')
  .reduceRegion({reducer: ee.Reducer.mean(), geometry: studyArea, scale: 30, maxPixels: 1e9}));
print('2023 NDVI mean:', comp2023.select('NDVI')
  .reduceRegion({reducer: ee.Reducer.mean(), geometry: studyArea, scale: 30, maxPixels: 1e9}));

// Tighter range visualization
Map.addLayer(ndviChange, {
  min: -0.15, max: 0.15,
  palette: ['red','orange','white','lightgreen','darkgreen']
}, 'NDVI Change (tight range)', true);
