// ================================================
// PROJECT: Deforestation Monitoring - Rangamati
// Script 01: Study Area + Data Exploration
// ================================================

// --- Study Area (Rangamati District) ---
var bgd = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
);
var rangamati = bgd.filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
var studyArea = rangamati.geometry();

// --- Landsat 5 (2000) Composite ---
var l5_2000 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2000-02-01', '2000-03-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 20))
  .median()
  .clip(studyArea);

// --- Landsat 8 (2023) Composite ---
var l8_2023 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2023-02-01', '2023-03-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 20))
  .median()
  .clip(studyArea);

// --- Map Display ---
Map.centerObject(studyArea, 10);

// Rangamati boundary
Map.addLayer(rangamati, {color: 'red'}, 'Rangamati District');

// Landsat 5 - 2000 (True Color)
Map.addLayer(l5_2000,
  {bands: ['SR_B3','SR_B2','SR_B1'], min: 7000, max: 13000},
  'L5 True Color 2000');

// Landsat 8 - 2023 (True Color)
Map.addLayer(l8_2023,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 7000, max: 13000},
  'L8 True Color 2023');

// --- Console Info ---
print('Study Area: Rangamati Hill District');
print('Area (km²):', studyArea.area().divide(1e6));
print('L5 images (2000):', ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(studyArea).filterDate('2000-02-01','2000-03-31').size());
print('L8 images (2023):', ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(studyArea).filterDate('2023-02-01','2023-03-31').size());
