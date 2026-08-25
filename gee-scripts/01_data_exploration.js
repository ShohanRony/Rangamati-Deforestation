// ============================================================
// PROJECT: Deforestation Monitoring - Rangamati
// Script 01a: NDVI Data Exploration (Landsat 5 vs Landsat 8)
// Researcher: Shohinur Pervez Shohan, RMSTU
// Purpose: Compare NDVI vegetation cover between 2000 and 2023
//          using an approximate bounding rectangle for the district.
// ============================================================

// Rangamati district - approximate bounding rectangle
var rangamati = ee.Geometry.Rectangle([91.90, 22.05, 92.55, 23.75]);

// ---- Landsat 5: 2000 ----
var l5_2000 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
.filterBounds(rangamati)
.filterDate('2000-02-01', '2000-03-31')
.filter(ee.Filter.lt('CLOUD_COVER', 30));

// ---- Landsat 8: 2023 ----
var l8_2023 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
.filterBounds(rangamati)
.filterDate('2023-02-01', '2023-03-31')
.filter(ee.Filter.lt('CLOUD_COVER', 30));

print('L5 images (2000 Feb-Mar):', l5_2000.size());
print('L8 images (2023 Feb-Mar):', l8_2023.size());

// NDVI 2000 vs 2023
var ndvi_2000 = l5_2000.median()
.normalizedDifference(['SR_B4', 'SR_B3'])
.rename('NDVI_2000');

var ndvi_2023 = l8_2023.median()
.normalizedDifference(['SR_B5', 'SR_B4'])
.rename('NDVI_2023');

// Map center: Rangamati Sadar (22 37'N, 92 12'E)
Map.setCenter(92.12, 22.62, 10);

Map.addLayer(ndvi_2000,
             {min: 0.1, max: 0.8, palette: ['red','yellow','green']},
             'NDVI 2000 (L5)');

Map.addLayer(ndvi_2023,
             {min: 0.1, max: 0.8, palette: ['red','yellow','green']},
             'NDVI 2023 (L8)');

// Display district boundary
Map.addLayer(rangamati, {color: 'blue'}, 'Rangamati Boundary');

print('Complete!');
