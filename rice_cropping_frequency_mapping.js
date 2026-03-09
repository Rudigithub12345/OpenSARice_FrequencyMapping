//Fatchurrachman, F., Rudiyanto, R., Ginting, F.I., Che Soh, N., Mohd Shah, R., Goh Eng Giap, S.,  
//Setiawan, B.I., Powell, V., Hsu, M., Schiller, S., Davitt, A., Setia, R., & Minasny, B. (2026). 
//Rice extent and cropping frequency dataset across South Asia derived from Sentinel-1 and Sentinel-2 time series. 

//email: rudiyanto@umt.edu.my



// Define Region of Interest (ROI)

var geometry_roi = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[88.99855407732021, 25.00089529766422],
          [88.99718078630458, 23.998228704290934],
          [90.00002655046474, 23.99901281505026],
          [90.00071319597255, 25.000895297664243]]], null, false);


var roi = geometry_roi;
Map.centerObject(roi);

// Load and process ESA WorldCover dataset
var dataset_woco = ee.ImageCollection("ESA/WorldCover/v100").first().clip(roi);
var crop_woco = dataset_woco.eq(40).or(dataset_woco.eq(90)); // 40 = cropland, 90 = built-up areas
crop_woco = crop_woco.updateMask(crop_woco);
Map.addLayer(crop_woco, {min: 0, max: 1, palette: ['white', 'red']}, "Cropland wo-co", false);

var land = dataset_woco.neq(80).clip(roi); // Exclude water bodies (80)
Map.addLayer(land, {}, 'Land Cover', false);

// Visualize ROI and compute area
Map.addLayer(roi, {color: 'blue'}, 'ROI', false);
print('Area of ROI (ha)', roi.area({maxError: 1}).divide(1e4).round());

// Sentinel-2 Image Collection Setup
var startDate1 = ee.Date('2020-01-01');
var endDate1 = ee.Date('2021-12-31');
var S2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filter(ee.Filter.date(startDate1, endDate1))
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 100)
  .filterBounds(roi);

print('Sample Sentinel-2 Images', S2.limit(10));
print('Sentinel-2 Satellites', S2.aggregate_array('SPACECRAFT_NAME'));

// Functions to calculate NDVI and NDSI
var addNDVI = function(image) {
  return image.addBands(image.normalizedDifference(['B8', 'B4']).rename('ndvi'));
};

var addNDSI = function(image) {
  return image.addBands(image.normalizedDifference(['B3', 'B11']).rename('ndsi'));
};

// Map NDVI computation across the dataset
var datasetS2 = S2.map(addNDVI);
print('Dataset with NDVI', datasetS2);

// Temporal NDVI Quality Mosaic
var month = 1;
var months = endDate1.difference(startDate1, 'month').divide(month).toInt();
var sequence = ee.List.sequence(0, months);
print('Time Sequence', sequence);

var sequence_s2 = sequence.map(function(num) {
  num = ee.Number(num);
  var Start_interval = startDate1.advance(num.multiply(month), 'month');
  var End_interval = startDate1.advance(num.add(1).multiply(month), 'month');
  var subset = datasetS2.filterDate(Start_interval, End_interval);
  var date_img = Start_interval.format('yyyy-MM-dd');
  return subset.qualityMosaic('ndvi')
               .set('system:time_start', Start_interval)
               .set({'system:time_index': date_img});
});

var byMonthYear = ee.ImageCollection.fromImages(sequence_s2).filter(ee.Filter.date(startDate1, endDate1));
print('Monthly NDVI Images', byMonthYear);

// Map NDSI computation across the dataset
var datasetNDSI = byMonthYear.map(addNDSI);
print('Dataset with NDSI', datasetNDSI);

// Convert NDVI and NDSI collections to multiband images
var multibandNDVI = datasetNDSI.select('ndvi').toBands().clip(roi);
var multibandNDSI = datasetNDSI.select('ndsi').toBands().clip(roi);

// Handle system:time_index renaming
var monList = byMonthYear.aggregate_array('system:time_index');
multibandNDVI = multibandNDVI.rename(monList).clip(roi);
multibandNDSI = multibandNDSI.rename(monList).clip(roi);

print('Multiband NDVI', multibandNDVI);
print('Multiband NDSI', multibandNDSI);

// Visualize NDVI on the map
Map.addLayer(multibandNDVI, {min: 0.25, max: 0.65}, 'Monthly S-2 NDVI', false);

// Sentinel-1 VH Image Collection Setup
var sentinel1_vh = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) // Filter for VH polarization
  .select('VH') // Select the VH band
  .filter(ee.Filter.eq('instrumentMode', 'IW')) // Use Interferometric Wide swath mode
  .filter(ee.Filter.eq('resolution_meters', 10)) // Ensure 10-meter resolution
  .filter(ee.Filter.date(startDate1, endDate1)) // Filter by date range
  .filter(ee.Filter.bounds(roi)); // Filter by the ROI
print('Filtered Sentinel-1 VH Collection', sentinel1_vh);

// Calculate number of months in the time range
var month = 1; // Define interval duration (1 month)
var months = endDate1.difference(startDate1, 'month').divide(month).toInt(); // Total intervals
print('Number of Months', months);

// Generate a sequence for monthly intervals
var sequence = ee.List.sequence(0, months);
print('Time Sequence', sequence);

// Function to create monthly composites
var sequence_s1 = sequence.map(function(num) {
  num = ee.Number(num); // Convert sequence number to EE number
  var Start_interval = startDate1.advance(num.multiply(month), 'month'); // Start of interval
  var End_interval = startDate1.advance(num.add(1).multiply(month), 'month'); // End of interval
  
  // Filter Sentinel-1 data for the interval and calculate the median
  var subset = sentinel1_vh.filterDate(Start_interval, End_interval);
  var date_img = Start_interval.format('yyyy-MM-dd'); // Label for the time index
  
  // Return the median image with metadata for time
  return subset.median()
               .set('system:time_start', Start_interval)
               .set({'system:time_index': date_img});
});

// Create ImageCollection from monthly composites
var byMonthYearS1a = ee.ImageCollection.fromImages(sequence_s1)
                                       .filter(ee.Filter.date(startDate1, endDate1));
print('Monthly Median VH Collection', byMonthYearS1a);


// Debugging comb_byMonthYears1
print("comb_byMonthYears1", byMonthYearS1a); // Assuming 'byMonthYearS1a' is the input
var listOfImages_S1 = byMonthYearS1a.toList(byMonthYearS1a.size()); // Convert to list for processing
print("List of Images (S1)", listOfImages_S1);

// Filter images with only one band
var newList = listOfImages_S1.map(function(ele) {
  ele = ee.Image(ele); // Ensure the element is an image
  var bandCount = ele.bandNames().size(); // Count the number of bands
  return ee.Algorithms.If(bandCount.eq(1), ele, null); // Retain images with 1 band
}).removeAll([null]); // Remove null elements
print("Filtered List (Single-band Images)", newList);

// Create a new ImageCollection from the filtered list
var comb_byMonthYears1 = ee.ImageCollection(newList);
print("Filtered ImageCollection", comb_byMonthYears1);

// Extract 'system:time_index' property for renaming bands
var monList = comb_byMonthYears1.aggregate_array('system:time_index');
print("Time Index List", monList);

// Convert the ImageCollection into a multiband image and rename bands
var multibands1 = comb_byMonthYears1.toBands().rename(monList).clip(roi);

// Add the multiband image to the map for visualization
Map.addLayer(multibands1, {min: -25, max: -15}, 'Monthly S-1 VH', true);



// Combine layers (NDVI, NDSI, VH) and mask using 'land' and 'crop_woco'
var combinedband = multibandNDVI
  .addBands(multibandNDSI)
  .addBands(multibands1)
  .updateMask(land)
  .updateMask(crop_woco);
print('Combined Band Image', combinedband);

// Extract band names from the combined image
var listBandNames = combinedband.bandNames();
print('Band Names in Combined Image', listBandNames);

// Sample points for clustering
var training = combinedband.sample({
  region: roi,
  scale: 10,
  tileScale: 4,
  numPixels: 3000,
  geometries: true
});
//Map.addLayer(training, {color: 'blue'}, 'Training Points', false);
print('Training Dataset', training);

// Train K-Means clusterer
var clusterer = ee.Clusterer.wekaKMeans(25).train({
  features: training,
  inputProperties: listBandNames
});

// Apply clustering to combined bands
var result_cluster = combinedband.select(listBandNames).cluster(clusterer).byte();
print('Clustered Image', result_cluster);

// Define remapping values
var clusters = ee.List.sequence(0, 29); // Original cluster IDs (0–29)
var values0 = ee.List.sequence(1, 30);  // Remapped values for the first reclassification
var values1 =   [0, 1, 0, 1,1,
                 1, 1, 1, 1,0,//10
                 0, 1, 1, 1,1,
                 1, 1, 1, 1,1,//20
                 0, 1, 1, 1,1,
                 0, 0, 0, 0,0];//30
                 
                 
var values2 =   [0, 2, 0, 1,2,
                 2, 2, 2, 1,0,//10
                 0, 2, 2, 1,2,
                 2, 1, 1, 2,2,//20
                 0, 2, 2, 2,1,
                 0, 0, 0, 0,0];//30                 

// First remapping (arbitrary values)
var remapped_cluster = result_cluster.remap(clusters, values0).byte();

// Second remapping (binary classification)
var remapped_cluster1 = remapped_cluster.remap(values0, values1).byte()
 // .rename('classification');
 // pattern
var remapped_cluster2 = remapped_cluster.remap(values0, values2).byte() 


// Visualization layers
Map.addLayer(remapped_cluster.randomVisualizer(), {}, 'Reclassified Clusters (Arbitrary)', true);
//Map.addLayer(remapped_cluster1.randomVisualizer(), {}, 'Reclassified Clusters (Binary)', true);

// Binary classification visualization
var binary_paddy = remapped_cluster1.updateMask(remapped_cluster1);
var pattern_paddy = remapped_cluster2.updateMask(remapped_cluster1);

Map.addLayer(binary_paddy, {min: 0, max: 1, palette: ['blue', 'red']}, 'Binary Classification', true);
Map.addLayer(pattern_paddy.randomVisualizer(), {}, 'Pattern Paddy', true);



// Combine NDVI multiband image with remapped clusters
var comb_ndvi_cluster = multibandNDVI.addBands(remapped_cluster);
print("Combined NDVI and Cluster", comb_ndvi_cluster);

// Define chart options for NDVI
var options_ndvi = {
  lineWidth: 1,
  pointSize: 2,
  hAxis: {title: 'Year-Month'},
  vAxis: {title: 'NDVI'},
  title: 'Sentinel-2 NDVI Spectra in Classified Regions'
};

// Create NDVI chart by class
// Make the chart, set the options.
var chart_class_ndvi = ui.Chart.image.byClass(
    comb_ndvi_cluster, 'remapped', roi, ee.Reducer.median(), 500)//, classNames, wavelengths)
    .setOptions(options_ndvi)
    .setChartType('ScatterChart');
//print("NDVI Spectra Chart", chart_class_ndvi);


// Combine NDSI multiband image with remapped clusters
var comb_ndsi_cluster = multibandNDSI.addBands(remapped_cluster);
print("Combined NDSI and Cluster", comb_ndsi_cluster);

// Define chart options for NDSI
var options_ndsi = {
  lineWidth: 1,
  pointSize: 2,
  hAxis: {title: 'Year-Month'},
  vAxis: {title: 'NDSI'},
  title: 'Sentinel-2 NDSI Spectra in Classified Regions'
};

// Create NDSI chart by class
var chart_class_ndsi = ui.Chart.image.byClass(
    comb_ndsi_cluster, 'remapped', roi, ee.Reducer.median(), 500) //, classNames, wavelengths)
    .setOptions(options_ndsi)
    .setChartType('ScatterChart');

// Print the NDSI chart
//print("NDSI Spectra Chart", chart_class_ndsi);





// Combine VH multiband image with remapped clusters
var comb_vh_cluster = multibands1.addBands(remapped_cluster);
print("Combined VH and Cluster", comb_vh_cluster);

// Define chart options for VH
var options_vh = {
  lineWidth: 1,
  pointSize: 2,
  hAxis: {title: 'Year-Month'},
  vAxis: {title: 'VH'},
  title: 'Sentinel-1 VH Spectra in Classified Regions'
};

// Create VH chart by class
// Make the chart, set the options.
var chart_class_vh = ui.Chart.image.byClass(
    comb_vh_cluster, 'remapped', roi, ee.Reducer.median(), 500)//, classNames, wavelengths)
    .setOptions(options_vh)
    .setChartType('ScatterChart');
//print("VH Spectra Chart", chart_class_vh);
//

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
    fontSize: '14px',
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
    fontSize: '16px',
    color: 'black'
  }
);

// Combine affiliation into a panel
var affiliationPanel = ui.Panel([affiliationLabel, affiliation]);

// Add components to the main panel
panel.add(label1);
panel.add(label2);
//panel.add(affiliationPanel);
panel.add(chart_class_ndvi);
panel.add(chart_class_ndsi);
panel.add(chart_class_vh);




