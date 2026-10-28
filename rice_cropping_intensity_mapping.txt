//Pribadi, A., Rudiyanto, Fatchurrachman, Ginting, F. I., Che Soh, N., Mohd Shah, R., Kok, S. C., 
//Goh, E. G., Setiawan, B. I., de Vos, K., Andriamananjara, A., Massawe, B. H. J., Supriatna, Powell, 
//V., Hsu, M., Schiller, S., Davitt, A., & Minasny, B. (2026). 
//A high resolution dataset of paddy rice distribution and cropping intensity across East Africa. 

//email: rudiyanto@umt.edu.my



// =============================================================================
// STEP 1: max-NDVI (qualityMosaic) monthly composites, three series whose start
// dates are shifted by 0 / 10 / 20 days, then merged = 30 windows.
// No cloud mask, no gap filling.
// Order: (1) NDVI per image -> (2) qualityMosaic on NDVI -> (3) MNDWI computed
// from the composite's own B3/B11 (i.e. from the max-NDVI pixel).
// Requires import: geometry_roi
// =============================================================================

// ------------------------------- CONFIG --------------------------------------
var START_DATE = ee.Date('2024-01-01'); // <- known start date (edit this)
var N_MONTHS = 12;                      // months per series
var OFFSETS  = [0, 10, 20];             // start shift (days) of each series
// 11 months x 3 offsets = 33 windows (use N_MONTHS = 10 for 30 windows)
 
// End derived: the LAST window (last month, last offset) ends here (exclusive)
var END_DATE = START_DATE.advance(OFFSETS[OFFSETS.length - 1], 'day')
                         .advance(N_MONTHS, 'month');

var roi = geometry_roi;
Map.centerObject(roi,10);
Map.setOptions('HYBRID');


var dataset_woco = ee.ImageCollection("ESA/WorldCover/v100").first();//.clip(roi);

var visualization = {
  bands: ['Map'],
};


Map.addLayer(dataset_woco, visualization, 'Landcover',0);


var crop_woco=dataset_woco.eq(40).or(dataset_woco.eq(30)).or(dataset_woco.eq(20));//.or(dataset_woco.eq(90)).or(dataset_woco.eq(95)).or(dataset_woco.eq(20)).or(dataset_woco.eq(30));//.or(dataset_woco.eq(10));//.or(dataset_woco.eq(30)).or(dataset_woco.eq(80));//.or(dataset_woco.eq(80))//40 cropland 50 builtup
var crop_woco=crop_woco.updateMask(crop_woco)
Map.addLayer(crop_woco, {min:0,max:1,palette:['white','red']}, "Cropland wo-co",0);

var land=dataset_woco.neq(80).clip(roi)

print('Period', START_DATE.format('YYYY-MM-dd'), END_DATE.format('YYYY-MM-dd'));

// ------------------------------- SENTINEL-2 ----------------------------------
// Only NDVI is computed before the composite; B3 and B11 are kept for MNDWI later
function addNDVI(img) {
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('ndvi');
  return img.select(['B3', 'B11']).toFloat().addBands(ndvi.toFloat());
}

var s2All = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filterBounds(roi)
  .filterDate(START_DATE, END_DATE)
  .map(addNDVI);

// Masked blank keeps band names when a window has no scenes
var blank = ee.Image.constant([0, 0, 0]).rename(['B3', 'B11', 'ndvi']).toFloat().selfMask();

// One function replaces the three copy-pasted blocks
function monthlySeries(offsetDays) {
  var s0 = START_DATE.advance(offsetDays, 'day');
  return ee.ImageCollection.fromImages(
    ee.List.sequence(0, N_MONTHS - 1).map(function (m) {
      m = ee.Number(m);
      var s = s0.advance(m, 'month');
      var e = s0.advance(m.add(1), 'month');
      return s2All.filterDate(s, e)
        .merge(ee.ImageCollection(blank))
        .qualityMosaic('ndvi')                       // step 2: max NDVI pixel
        .set({
          date: s.format('yyyy-MM-dd'),
          'system:time_start': s.millis()
        });
    }));
}

// Merge the three series and sort by time
var composites = monthlySeries(OFFSETS[0]);
for (var k = 1; k < OFFSETS.length; k++) {
  composites = composites.merge(monthlySeries(OFFSETS[k]));
}
composites = composites.sort('system:time_start');

// Step 3: MNDWI applied AFTER qualityMosaic, on the composite's B3 / B11
composites = composites.map(function (img) {
  var mndwi = img.normalizedDifference(['B3', 'B11']).rename('mndwi');
  return img.addBands(mndwi);
});

print('composites (30)', composites);
print('window start dates', composites.aggregate_array('date'));

// ------------------------------- OUTPUT --------------------------------------
var dates = composites.aggregate_array('date');
var ndviNames = dates.map(function (d) { return ee.String('ndvi_').cat(d); });

var ndviStack = composites.select('ndvi').toBands().rename(ndviNames)
  .updateMask(land).clip(roi);
print('ndviStack (30 bands)', ndviStack);

var mndwiNames = dates.map(function (d) { return ee.String('mndwi_').cat(d); });
var mndwiStack = composites.select('mndwi').toBands().rename(mndwiNames)
  .updateMask(land).clip(roi);
print('mndwiStack (30 bands)', mndwiStack);

var ndviVis = {min: 0.2, max: 0.8, palette: ['red', 'orange', 'yellow', 'green']};
Map.addLayer(ndviStack.reduce(ee.Reducer.max()), ndviVis, 'max NDVI (all windows)', 1);
Map.addLayer(mndwiStack.reduce(ee.Reducer.mean()),
  {min: -0.6, max: 0.3, palette: ['brown', 'white', 'blue']}, 'mean MNDWI (all windows)', 0);
Map.addLayer(ee.Image(composites.toList(30).get(0)).select('ndvi').clip(roi),
  ndviVis, 'NDVI first window', 0);
  
  
///sentinel-1
// =============================================================================
// STEP 1: max-NDVI (qualityMosaic) monthly composites, three series whose start
// dates are shifted by 0 / 10 / 20 days, then merged (N_MONTHS x 3 windows).
// No cloud mask, no gap filling, no S1, no clustering yet.
// Order: (1) NDVI per image -> (2) qualityMosaic on NDVI -> (3) MNDWI computed
// from the composite's own B3/B11 (i.e. from the max-NDVI pixel).
// Then: Sentinel-1 median VH on the SAME windows (section at the end).
// Requires import: geometry_roi
// =============================================================================
/*
// ------------------------------- CONFIG --------------------------------------
var START_DATE = ee.Date('2024-01-01'); // <- known start date (edit this)
var N_MONTHS = 12;                      // months per series
var OFFSETS  = [0, 10, 20];             // start shift (days) of each series
// 11 months x 3 offsets = 33 windows (use N_MONTHS = 10 for 30 windows)

// End derived: the LAST window (last month, last offset) ends here (exclusive)
var END_DATE = START_DATE.advance(OFFSETS[OFFSETS.length - 1], 'day')
                         .advance(N_MONTHS, 'month');

var roi = geometry_roi;
//Map.centerObject(roi, 8);
Map.setOptions('HYBRID');
*/
var land = ee.ImageCollection('ESA/WorldCover/v100').first().neq(80);

print('Period', START_DATE.format('YYYY-MM-dd'), END_DATE.format('YYYY-MM-dd'));

// ------------------------------- SENTINEL-2 ----------------------------------
// Only NDVI is computed before the composite; B3 and B11 are kept for MNDWI later
function addNDVI(img) {
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('ndvi');
  return img.select(['B3', 'B11']).toFloat().addBands(ndvi.toFloat());
}

var s2All = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filterBounds(roi)
  .filterDate(START_DATE, END_DATE)
  .map(addNDVI);

// Masked blank keeps band names when a window has no scenes
var blank = ee.Image.constant([0, 0, 0]).rename(['B3', 'B11', 'ndvi']).toFloat().selfMask();

// One function replaces the three copy-pasted blocks
function monthlySeries(offsetDays) {
  var s0 = START_DATE.advance(offsetDays, 'day');
  return ee.ImageCollection.fromImages(
    ee.List.sequence(0, N_MONTHS - 1).map(function (m) {
      m = ee.Number(m);
      var s = s0.advance(m, 'month');
      var e = s0.advance(m.add(1), 'month');
      return s2All.filterDate(s, e)
        .merge(ee.ImageCollection(blank))
        .qualityMosaic('ndvi')                       // step 2: max NDVI pixel
        .set({
          date: s.format('yyyy-MM-dd'),
          'system:time_start': s.millis()
        });
    }));
}

// Merge the three series and sort by time
var composites = monthlySeries(OFFSETS[0]);
for (var k = 1; k < OFFSETS.length; k++) {
  composites = composites.merge(monthlySeries(OFFSETS[k]));
}
composites = composites.sort('system:time_start');

// Step 3: MNDWI applied AFTER qualityMosaic, on the composite's B3 / B11
composites = composites.map(function (img) {
  var mndwi = img.normalizedDifference(['B3', 'B11']).rename('mndwi');
  return img.addBands(mndwi);
});

print('composites', composites);
print('window start dates', composites.aggregate_array('date'));

// ------------------------------- OUTPUT --------------------------------------
var dates = composites.aggregate_array('date');
var ndviNames = dates.map(function (d) { return ee.String('ndvi_').cat(d); });

var ndviStack = composites.select('ndvi').toBands().rename(ndviNames)
  .updateMask(land).clip(roi);
print('ndviStack', ndviStack);

var mndwiNames = dates.map(function (d) { return ee.String('mndwi_').cat(d); });
var mndwiStack = composites.select('mndwi').toBands().rename(mndwiNames)
  .updateMask(land).clip(roi);
print('mndwiStack', mndwiStack);

var ndviVis = {min: 0.2, max: 0.8, palette: ['red', 'orange', 'yellow', 'green']};
Map.addLayer(ndviStack.reduce(ee.Reducer.max()), ndviVis, 'max NDVI (all windows)', 1);
Map.addLayer(mndwiStack.reduce(ee.Reducer.mean()),
  {min: -0.6, max: 0.3, palette: ['brown', 'white', 'blue']}, 'mean MNDWI (all windows)', 0);
Map.addLayer(ee.Image(composites.first()).select('ndvi').updateMask(land).clip(roi),
  ndviVis, 'NDVI first window', 0);


// =============================================================================
// SENTINEL-1: median VH on the same windows (same START_DATE, N_MONTHS, OFFSETS)
// =============================================================================
var s1All = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('resolution_meters', 10))
  // .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) // optional: single pass
  .filterBounds(roi)
  .filterDate(START_DATE, END_DATE)
  .select('VH');                                   // already in dB

// Masked blank keeps the VH band (and the window) when no S1 scene exists.
// This replaces the old bandNames().size() / removeAll([0]) workaround, which
// DROPPED empty windows and broke alignment with the NDVI/MNDWI windows.
var blankVH = ee.Image.constant(0).rename('VH').toFloat().selfMask();

function s1MonthlySeries(offsetDays) {
  var s0 = START_DATE.advance(offsetDays, 'day');
  return ee.ImageCollection.fromImages(
    ee.List.sequence(0, N_MONTHS - 1).map(function (m) {
      m = ee.Number(m);
      var s = s0.advance(m, 'month');
      var e = s0.advance(m.add(1), 'month');
      return s1All.filterDate(s, e)
        .merge(ee.ImageCollection(blankVH))
        .median()
        .rename('VH')
        .set({
          date: s.format('yyyy-MM-dd'),
          'system:time_start': s.millis()
        });
    }));
}

var s1Windows = s1MonthlySeries(OFFSETS[0]);
for (var j = 1; j < OFFSETS.length; j++) {
  s1Windows = s1Windows.merge(s1MonthlySeries(OFFSETS[j]));
}
s1Windows = s1Windows.sort('system:time_start');

print('s1Windows', s1Windows);
print('S1 window dates (must equal S2 dates)', s1Windows.aggregate_array('date'));

var vhNames = s1Windows.aggregate_array('date')
  .map(function (d) { return ee.String('vh_').cat(d); });

var vhStack = s1Windows.select('VH').toBands().rename(vhNames)
  .updateMask(land).clip(roi);

// Normalise dB to 0-1 (same -25..-10 range as the original), clamped
var VH_MIN = -25, VH_MAX = -10;
var vhNorm = vhStack.subtract(VH_MIN).divide(VH_MAX - VH_MIN).clamp(0, 1);

print('vhStack', vhStack);
print('vhNorm', vhNorm);
Map.addLayer(vhNorm.reduce(ee.Reducer.min()).addBands(vhNorm.reduce(ee.Reducer.median()).addBands(vhNorm.reduce(ee.Reducer.max()))), {min: 0.2, max: 0.6}, 'norm VH (the first three windows)', 0);
Map.addLayer(vhStack.reduce(ee.Reducer.median()), {min: -25, max: -12}, 'median VH (all windows)', 0);
Map.addLayer(vhStack.select(0), {min: -25, max: -12}, 'VH first window', 0);

var combinedband=ndviStack.addBands(mndwiStack).addBands(vhStack).updateMask(crop_woco).updateMask(land)
//var combinedband=ndviStack.addBands(mndwiStack).addBands(vhNorm).updateMask(crop_woco).updateMask(land)

print('comb_img',combinedband)

var listBandNames=combinedband.bandNames();
print('listBandNames',listBandNames);


var N_CLUSTERS = 15;
var N_SAMPLES  = 3000;
var SCALE      = 10;
var SEED       = 42;
var MAX_GAP    = 3;     // max gap (in windows) to interpolate; windows are ~10 d apart
// After inspecting the clusters, list the IDs to keep as class 1, e.g. [2, 7]
var TARGET_CLUSTERS = [];

// ------------------------------- K-MEANS -------------------------------------
var training = combinedband.sample({
  region: roi,
  scale: SCALE,
  numPixels: N_SAMPLES,
  seed: SEED,
  tileScale: 8,
  geometries:true
});
print('training size', training.size());
 
var clusterer = ee.Clusterer.wekaKMeans({
  nClusters: N_CLUSTERS,
  init: 1,               // k-means++
  maxIterations: 100,
  seed: SEED
}).train({
  features: training,
  inputProperties: listBandNames
});
 
var result_cluster = combinedband.select(listBandNames).cluster(clusterer)
    .add(1)                 // clusters 1..N_CLUSTERS instead of 0..N_CLUSTERS-1
  .rename('cluster').byte();
Map.addLayer(result_cluster.randomVisualizer(), {}, 'k-means clusters', 1);
 
print('clusterer',clusterer)

var result_cluster_point = training.select(listBandNames).cluster(clusterer)
  .map(function (f) {
    return f.set('cluster', ee.Number(f.get('cluster')).add(1));
  });
print('result_cluster_point', result_cluster_point);
print('result_cluster_point',result_cluster_point);

var clusterIds = ee.List.sequence(1, N_CLUSTERS);

var values1 =   [0, 0, 0, 2,0,//5
                 2, 0, 0, 2,2,//10
                 0, 0, 2, 2,0//15
                // 0, 0, 2, 0,0,//19
               //  2, 0, 0, 2,0,
               //  0, 0, 0, 0,0,//];//30
               //  0, 0, 0, 0,0,
              //   0, 0, 0, 0,0
              ];
                               


                                  
var remapped_cluster = result_cluster.remap(clusterIds, values1).byte().selfMask(); 
Map.addLayer(remapped_cluster, {min:0,max:2,palette:['white','blue','red']}, 'remapped_cluster_intensity',1);  

Export.image.toDrive({
  image:remapped_cluster.toByte(), 
  description:'ExportToDrive_rice_intensity', 
  folder:'GEE_EastAfricaRice', 
  fileNamePrefix:'RiceIntensity', 
  region:roi, 
  scale:10, 
  crs:'EPSG:4326', 
   maxPixels:1e13, 

  })

Export.image.toAsset({
  image:remapped_cluster.toByte(),  
  description:'ExportToAssets_rice_intensity', 
  assetId:'RiceIntensity', 
  region:roi, 
  scale:10, 
  crs:'EPSG:4326', 
  maxPixels:1e13})


// =============================================================================
// Median profile per cluster for NDVI, MNDWI and VH (one chart each)
// Paste AFTER:  var result_cluster_point = training.select(listBandNames).cluster(clusterer)
// Uses: ndviStack, mndwiStack, vhStack (band names end with the window date,
//       e.g. 'ndvi_2024-01-11'), and the 'cluster' property on the points.
// =============================================================================


// Create a panel next to the map with defined width
var panel = ui.Panel({style: {width: '550px'}});
ui.root.add(panel);

// Add a title label to the panel
var label1 = ui.Label({
  value: 'SENTINEL-2 & SENTINEL-1 Explorer: Paddy Clustering',
  style: {
    color: 'white',
    backgroundColor: 'blue',
    fontWeight: 'bold',
    fontSize: '16px',
    border: '1px solid black',
    padding: '5px 5px 5px 5px',
    margin: '12px 0px 0px 8px'
  }
});

// Add a developer label
var label2 = ui.Label({
  value: 'Developed by Rudiyanto',
  style: {
    fontSize: '16px',
    color: 'black',
    padding: '5px'
  }
});

 var label3 = ui.Label({
  value: '*Email: rudiyanto@umt.edu.my',
  style: {
    fontSize: '12px',
    color: 'black',
   // fontWeight: 'bold',
    padding: '2px',
  }
});

// Add affiliation and instructions
var affiliationLabel = ui.Label('Instructions', {
  fontWeight: 'bold',
  padding: '3px'
});
var affiliation = ui.Label(
  "Program of Crop Science\n" +
  "Faculty of Fisheries and Food Science\n" +
  "Universiti Malaysia Terengganu\n" +
  "Email: rudiyanto@umt.edu.my\n", 
  {
    whiteSpace: 'pre',
    padding: '5px',
    fontSize: '14px',
    color: 'black'
  }
);

// Combine affiliation into a panel
//var affiliationPanel = ui.Panel([affiliationLabel, affiliation]);
var affiliationPanel = ui.Panel([affiliation]);
// Add components to the main panel
panel.add(label1);
panel.add(label2);
panel.add(affiliationPanel);
//panel.add(label2);
//panel.add(chart_class_ndvi);
//panel.add(chart_class_ndsi);
//panel.add(chart_class_vh);










var SHOW_SD = false;   // true = also draw an SD chart for each index

var INDICES = [
  {stack: ndviStack,  name: 'NDVI',     vMin: 0,  vMax: 1},    // fixed axis
  {stack: mndwiStack, name: 'MNDWI',     vMin: -0.6,  vMax: 0.2},                           // auto axis
  {stack: vhStack,    name: 'VH (dB)',     vMin: -24,  vMax: -14}                          // auto axis
];

// x-axis = real dates (yyyy-MM-dd). A plain string axis makes the chart fail, so the
// table is a DataTable literal with a column of type 'date'.
function toDateLiteral(name) {
  // name ends with 'YYYY-MM-DD' -> 'Date(Y, M-1, D)' (month is 0-based)
  var p = name.slice(-10).split('-');
  return 'Date(' + Number(p[0]) + ',' + (Number(p[1]) - 1) + ',' + Number(p[2]) + ')';
}

function drawChart(names, g, key, title, vTitle, vMin, vMax) {
  var cols = [{id: 'date', label: 'window start', type: 'date'}];
  g.forEach(function (c) {
    cols.push({id: 'c' + Math.round(c.group), label: 'cluster ' + Math.round(c.group), type: 'number'});
  });

  var rows = [];
  for (var i = 0; i < names.length; i++) {
    var cells = [{v: toDateLiteral(names[i])}];
    g.forEach(function (c) { cells.push({v: c[key][i]}); });
    rows.push({c: cells});
  }

  var vAxis = {title: vTitle};
  if (vMin !== undefined) vAxis.viewWindow = {min: vMin, max: vMax};

 // print
  panel.add(ui.Chart({cols: cols, rows: rows}, 'LineChart', {
    title: title,
    hAxis: {title: 'window start date', format: 'yyyy-MM-dd',
            slantedText: true, slantedTextAngle: 60},
    vAxis: vAxis,
    lineWidth: 2,
    pointSize: 3,
    interpolateNulls: true
  }));
}

function clusterProfile(cfg) {
  var bands = cfg.stack.bandNames();
  var bandsWithClass = bands.add('cluster');

  // median + stdDev on the same inputs, one pair per band, grouped by cluster
  var reducer = ee.Reducer.median()
    .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
    .repeat(bands.length())
    .group(bandsWithClass.indexOf('cluster'));

  var stats = result_cluster_point.reduceColumns({
    selectors: bandsWithClass,
    reducer: reducer
  });
  var groups = ee.List(stats.get('groups'));

  bands.evaluate(function (names) {
    groups.evaluate(function (g) {
      drawChart(names, g, 'median', 'Median ' + cfg.name + ' per cluster',
                cfg.name, cfg.vMin, cfg.vMax);
      if (SHOW_SD) {
        drawChart(names, g, 'stdDev', 'SD of ' + cfg.name + ' per cluster',
                  'SD ' + cfg.name);
      }
    });
  });
}

INDICES.forEach(clusterProfile);

print('points per cluster', result_cluster_point.aggregate_histogram('cluster'));


