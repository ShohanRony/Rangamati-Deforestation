// ============================================================================
// Rangamati Land-Cover Change, 1993-2023
// Script 10 - Per-Epoch Landsat Data-Quality Check (PRODUCTION)
// ----------------------------------------------------------------------------
// STATUS: production. Diagnostic used to test (and rule out) input-image
//   quality as an explanation for the 2003-2008-2018 forest-composition
//   swing reported in the manuscript.
//
// Purpose:  Report scene count, mean CLOUD_COVER metadata, and post-QA-mask
//           valid-pixel coverage for each of the six epochs within the
//           1 Jan - 30 Apr compositing window, to confirm compositing
//           input quality is consistent across epochs.
// Inputs:   projects/crypto-hallway-405211/assets/BGD_adm2 (private asset,
//           see docs/REPRODUCIBILITY.md); Landsat 5/7/8 C02 T1_L2
//           collections.
// Outputs:  Console prints (scene count, mean cloud cover %, valid-pixel
//           fraction) per epoch.
// Depends:  none.
// ============================================================================

var studyArea = ee.FeatureCollection(
  'projects/crypto-hallway-405211/assets/BGD_adm2'
).filter(ee.Filter.eq('NAME_2', 'Parbattya Chattagram'));
var roi = studyArea.geometry();

var epochs = [
  {year: 1993, col: 'LANDSAT/LT05/C02/T1_L2'},
  {year: 1998, col: 'LANDSAT/LT05/C02/T1_L2'},
  {year: 2003, col: 'LANDSAT/LE07/C02/T1_L2'},
  {year: 2008, col: 'LANDSAT/LT05/C02/T1_L2'},
  {year: 2018, col: 'LANDSAT/LC08/C02/T1_L2'},
  {year: 2023, col: 'LANDSAT/LC08/C02/T1_L2'},
];

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

epochs.forEach(function(ep) {
  var rawCol = ee.ImageCollection(ep.col)
    .filterBounds(roi)
    .filterDate(ep.year + '-01-01', ep.year + '-04-30');

  var sceneCount = rawCol.size();

  var maskedCol = rawCol.map(maskAndScale);

  // valid-pixel fraction of the median composite over the ROI
  var validMask = maskedCol.select('SR_B1').count(); // count of unmasked obs per pixel
  var validFrac = validMask.gt(0).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 100,
    maxPixels: 1e13,
    tileScale: 4
  });

  // mean CLOUD_COVER metadata across scenes (informational only, not used as a filter)
  var meanCloud = rawCol.aggregate_mean('CLOUD_COVER');

  print('Epoch ' + ep.year + ':',
    'scenes=', sceneCount,
    'meanCloudCoverPct=', meanCloud,
    'validPixelFraction=', validFrac.get('SR_B1'));
});
