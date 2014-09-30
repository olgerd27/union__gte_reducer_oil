//  PERFORMING A DEFINITION OF GTE's and OIL's SYSTEMS CHARACTERISTICS
tic;
xdel(winsid());
clear;
stacksize(5e7);
warning('off');
printf("*********************\n");
printf("* START application *\n");
printf("*********************\n");

// Additional definitions for start programm execution
TRUE = %T;      FALSE = %F;   // true-false definition
GTE_OIL = 1;    RED_OIL = 2;  // 1 - GTE diagnostics oil's system, 2 - reducer diagnostics oil's system

function path = getExecScriptPath()
//***************************************************************
// Function for obtain the current execution script file path   *
//***************************************************************
  [u, t, n] = file();
  index = grep(n, "/(?:.*\.sci|.*\.sce)$/", "r");
  path = fileparts(n(index(1)));
endfunction


//  INITIAL DATA
diag_sys = RED_OIL;       // PKSTD diagnostics system
gte_numb = 4;             // The GTE DA91 number

sectorLength = 1800;       // Length of splitting sectors
sectorShift = 300;        // Shift of sector with length sectorLength in every main cycle iteration
modelLength = 300;        // Length of the data model for performing of the forecasting temperature parameters on steady modes
forecastInterval = 300;   // Interval for the forecasting temperature parameters on steady modes

importSteady = TRUE;      // Import calculated steady mode points to the text file: TRUE - perform, FALSE - don't perform
if ~importSteady
  if diag_sys == GTE_OIL
    UGt_strange = 10;     // Settings Gt strange for defining the steady modes of GTE's work
    Un2_xx = 5500;        // Settings XX by n2 parameter for defining the steady modes of GTE's work
  else
    Ungv_strange = 5;     // Settings ngv strange for defining the steady modes of reducer's work
    Ungv_min = 40;        // Settings ngv min value for defining the steady modes of reducer's work
  end

  Udtm_valid_min = 0;     // Settings dtm min value for defining the invalid points
  Udtm_valid_max = 80;    // Settings dtm max value for defining the invalid points
end

if diag_sys == RED_OIL
  Nnom = 20020;           // The reducer power on the nominal reducer's work mode
  ngv_nom = 240;          // Rotation speed of the reducer outlet shaft on the nominal reducer's work mode
end

plotGraphs = TRUE;          // Graphs plot: TRUE - perform, FALSE - don't perform
plotGraphsSameWin = TRUE;   // Plot graphs, some quantity of that is placed on a same window: TRUE - perform, FALSE - don't perform

exportSteadyPoints  = FALSE; // Export steady mode points in the text file: TRUE - perform, FALSE - don't perform
exportResToTxtFile  = TRUE; // Export the oil's characteristics points values in a text file:  TRUE - perform, FALSE - don't perform
exportResToImgFiles = TRUE; // Export plotted the oil's characters points values in a graphics files with "png" extension:
                            // TRUE - perform, FALSE - don't perform

// The initial characteristics Ngte = f(p2) from a set of thermodynamics characteristics
p2_init = [1.795; 2.251; 2.702; 3.380; 4.147; 5.111; 6.322; 7.049; 7.775; 8.408;
           9.152; 9.836; 10.601; 12.189; 13.692; 15.132; 16.397; 17.560; 18.671; 19.782];
Ngte_init = [110.4; 273.4; 434.4; 676.5; 950.2; 1294.5; 2052.5; 2643.4; 3234.2; 3749.1;
             4353.8; 5098.5; 5958.6; 7823; 9869.7; 11981; 14013.4; 16015.3; 18017.2; 20019.1];

// Indexes of the oil's parameters
INIT_VALUE = 0;  index_in = INIT_VALUE; index_out = INIT_VALUE; // initialization
if diag_sys == GTE_OIL
  index_in = 1;  index_per = 2;  index_tkvd = 3;  index_tnd = 4;  index_tv = 5;  index_out = 6;
else
  index_in = 1;      index_out = 2;
  index_tz01a =  3;  index_tz01b =  4;  index_tz02a =  5;  index_tz02b =  6;  index_tz02c =  7;  index_tz02d =  8;
  index_tz03a =  9;  index_tz03b = 10;  index_tz03c = 11;  index_tz03d = 12;  index_tz04a = 13;  index_tz04b = 14;
  index_tz05a = 15;  index_tz05b = 16;  index_tz06a = 17;  index_tz06b = 18;  index_tz07a = 19;  index_tz07b = 20;
  index_tz08a = 21;  index_tz08b = 22;  index_tz09a = 23;  index_tz09b = 24;  index_tz10a = 25;  index_tz10b = 26;
  index_tz11a = 27;  index_tz11b = 28;
end

// Polynomial powers that describe the oil's characteristics dtX = f(Ngte)
if diag_sys == GTE_OIL
  polynPow(index_per - index_in) = 2;   // dtm_per
  polynPow(index_tkvd - index_in) = 2;  // dtm_tkvd
  polynPow(index_tnd - index_in) = 2;   // dtm_tnd
  polynPow(index_tv - index_in) = 2;    // dtm_tv
  polynPow(index_out - index_in) = 2;   // dtm_gte_out
else
  polynPow(index_out - index_in) = 1;   // dtm_red_out
  polynPow(index_tz01a - index_in) = 1; // dtz1a
  polynPow(index_tz01b - index_in) = 1; // dtz1b
  polynPow(index_tz02a - index_in) = 1; // dtz2a
  polynPow(index_tz02b - index_in) = 1; // dtz2b
  polynPow(index_tz02c - index_in) = 1; // dtz2c
  polynPow(index_tz02d - index_in) = 1; // dtz2d
  polynPow(index_tz03a - index_in) = 1; // dtz3a
  polynPow(index_tz03b - index_in) = 1; // dtz3b
  polynPow(index_tz03c - index_in) = 1; // dtz3c
  polynPow(index_tz03d - index_in) = 1; // dtz3d
  polynPow(index_tz04a - index_in) = 1; // dtz4a
  polynPow(index_tz04b - index_in) = 1; // dtz4b
  polynPow(index_tz05a - index_in) = 1; // dtz5a
  polynPow(index_tz05b - index_in) = 1; // dtz5b
  polynPow(index_tz06a - index_in) = 1; // dtz6a
  polynPow(index_tz06b - index_in) = 1; // dtz6b
  polynPow(index_tz07a - index_in) = 1; // dtz7a
  polynPow(index_tz07b - index_in) = 1; // dtz7b
  polynPow(index_tz08a - index_in) = 1; // dtz8a
  polynPow(index_tz08b - index_in) = 1; // dtz8b
  polynPow(index_tz09a - index_in) = 1; // dtz9a
  polynPow(index_tz09b - index_in) = 1; // dtz9b
  polynPow(index_tz10a - index_in) = 1; // dtz10a
  polynPow(index_tz10b - index_in) = 1; // dtz10b
  polynPow(index_tz11a - index_in) = 1; // dtz11a
  polynPow(index_tz11b - index_in) = 1; // dtz11b
end

if ~importSteady
  // Setting archives path, names and extension
  filesArchive = [// 1 etap PI
                  'mo_2013_1_3_23_0_0';   'mo_2013_1_3_9_52_48';  'mo_2013_1_4_9_34_5'; 'mo_2013_1_5_9_27_47';
                  'mo_2013_1_8_11_10_49'; 'mo_2013_1_8_16_25_21';
                  // 2 etap PI + PSI
                  'mo_2013_2_14_10_1_0';  'mo_2013_2_25_13_24_20';  'mo_2013_2_25_18_44_48'; 'mo_2013_2_25_9_58_3';
                  'mo_2013_2_26_9_28_29'; 'mo_2013_3_11_12_48_00'];

  //path_archives = "/media/oleg/users/Oleg/work_zm/export/GTA_M56/Archivs/GTE_DA91_#4_5/GTE_DA91_#4/all";
  path_archives = "D:\work\GTA_M56\Archivs\GTE_DA91_#4_5\GTE_DA91_#4\all";
end

// Names of oil's temperatures
if diag_sys == GTE_OIL
  t_names(index_in) = 'tm_gte_in';  t_names(index_per) = 'tm_per';
  t_names(index_tkvd) = 'tm_tkvd';  t_names(index_tnd) = 'tm_tnd';
  t_names(index_tv) = 'tm_tv';      t_names(index_out) = 'tm_gte_out';
else
  t_names(index_in) = 'tm_red_in';  t_names(index_out) = 'tm_red_out';
  t_names(index_tz01a) = 'tz01a';   t_names(index_tz01b) = 'tz01b';
  t_names(index_tz02a) = 'tz02a';   t_names(index_tz02b) = 'tz02b';
  t_names(index_tz02c) = 'tz02c';   t_names(index_tz02d) = 'tz02d';
  t_names(index_tz03a) = 'tz03a';   t_names(index_tz03b) = 'tz03b';
  t_names(index_tz03c) = 'tz03c';   t_names(index_tz03d) = 'tz03d';
  t_names(index_tz04a) = 'tz04a';   t_names(index_tz04b) = 'tz04b';
  t_names(index_tz05a) = 'tz05a';   t_names(index_tz05b) = 'tz05b';
  t_names(index_tz06a) = 'tz06a';   t_names(index_tz06b) = 'tz06b';
  t_names(index_tz07a) = 'tz07a';   t_names(index_tz07b) = 'tz07b';
  t_names(index_tz08a) = 'tz08a';   t_names(index_tz08b) = 'tz08b';
  t_names(index_tz09a) = 'tz09a';   t_names(index_tz09b) = 'tz09b';
  t_names(index_tz10a) = 'tz10a';   t_names(index_tz10b) = 'tz10b';
  t_names(index_tz11a) = 'tz11a';   t_names(index_tz11b) = 'tz11b';
end
//=============================================================================================================================

// LOADING additional files with functions
path_sourceFiles = getExecScriptPath();
names_sourceFiles = [
                    "add_functions.sci";
                    "in_out_functions.sci";
                    "special_functions.sci";
                    ];
for i = 1 : size(names_sourceFiles, 'r')
  exec(path_sourceFiles + names_sourceFiles(i)); // loading functionality from the script files
end
//=============================================================================================================================

// INITIALIZATION
// index_in and index_out
if (index_in == INIT_VALUE) | (index_out == INIT_VALUE)
  printf("[INFO]: The ""index_in"" or ""index_out"" is not defined in the INITIAL DATA section\n");
  return;
end

// params indexes arrays 
if diag_sys == GTE_OIL
  params_indexes = [index_in; index_per; index_tkvd; index_tnd; index_tv; index_out];
else
  params_indexes = [index_in; index_out; 
    index_tz01a; index_tz01b; index_tz02a; index_tz02b; index_tz02c; index_tz02d; index_tz03a; index_tz03b; 
    index_tz03c; index_tz03d; index_tz04a; index_tz04b; index_tz05a; index_tz05b; index_tz06a; index_tz06b;
    index_tz07a; index_tz07b; index_tz08a; index_tz08b; index_tz09a; index_tz09b; index_tz10a; index_tz10b;
    index_tz11a; index_tz11b];
end

count_tmParams = length(params_indexes); // quantity of the temperature
count_dtmParams = count_tmParams - 1; // quantity of the temperatures delta's

for i = 1 : count_dtmParams
  dt_names(i) = 'd' + t_names(i + 1); // names of the temperatures delta's
end

// Structure for storing parameters data, readed from an archive file
if diag_sys == GTE_OIL
  index_t0 = 1;  index_Gt = 2;  index_reg = 3;  index_n2 = 4;  index_tm = 5; // parameters indexes
  params(index_t0) = struct('name', 't0', 'archIndexStart', 16, 'archIndexEnd', 18, 'data', []);
  params(index_Gt) = struct('name', 'Gt', 'archIndexStart', 120, 'archIndexEnd', 120, 'data', []);
  params(index_reg) = struct('name', 'p2', 'archIndexStart', 38, 'archIndexEnd', 38, 'data', []);
  params(index_n2) = struct('name', 'n2', 'archIndexStart', 10, 'archIndexEnd', 10, 'data', []);
  params(index_tm) = struct('name', t_names, 'archIndexStart', 47, 'archIndexEnd', 52, 'data', []);
else
  index_reg = 1;  index_tm = 2; // parameters indexes
  params(index_reg) = struct('name', 'ngv', 'archIndexStart', 14, 'archIndexEnd', 14, 'data', []);
  params(index_tm) = struct('name', t_names, 'archIndexStart', 53, 'archIndexEnd', 80, 'data', []);
end

colors = [1, 2, 3, 5, 19, 16, 27, 22, 13, 6, 9, 32, 28, 21, 25, 23, 26, 17];

// INITIAL DATA FOR EXPORT
// Relative path's for data saving
path_dataRltv = "data"; // ralative path for all external data storing
path_intRltv = "int"; // relative path for saving internal data, that is need for programm work
path_resRltv = "out"; // relative path for saving results data
path_GTERltv = "gte"; // relative path for saving the GTE's results characteristics
path_reducerRltv = "reducer"; // relative path for saving the reducer's results characteristics
path_steadyRltv = "steady_modes"; // relative path for saving the steady mode points values
path_resTxtRltv = "txt"; // relative path for saving the text result
path_resImageRltv = "images"; // relative path for saving the images
// Define the result relative path in accordance with type of current PKSTD diagnostics system
if diag_sys == GTE_OIL
  path_diagSysRltv = path_GTERltv;
else
  path_diagSysRltv = path_reducerRltv;
end
sep = filesep(); // the dirs separator
// Results paths
// get a root path
indexes_sep = strindex(path_sourceFiles, sep);
index_last_ch = indexes_sep(length(indexes_sep) - 1) - 1; // index to last character before the last dir: [THIS INDEX]/[DIR NAME]/
path_root = part(path_sourceFiles, 1 : index_last_ch);
// forming the aim absolute paths
path_data = path_root + sep + path_dataRltv; // absolute path for all external data storing
path_int = path_data + sep + path_intRltv; // absolute path for saving the internal data
path_res = path_data + sep + path_resRltv; // absolute path for saving the results data
path_steady = path_int + sep + path_diagSysRltv + sep + path_steadyRltv; // the steady modes points absolute path
path_resTxt = path_res + sep + path_diagSysRltv + sep + path_resTxtRltv; // the text results absolute path
path_resImage = path_res + sep + path_diagSysRltv + sep + path_resImageRltv; // the graphic results absolute path
// Files extensions
ext_archive = 'txt'; // archives in-files extension
ext_steady = 'dat'; // steady out- and in-files extension
ext_out_images = 'png'; // images out-files extension
ext_out_txt = 'rez'; // text out-files extension
// Identification current calculation: [gte_numb]_[sectorLength]_[sectorShift]_[modelLength]_[forecastInterval]
str_currCalcIdentif = 'gn=' + string(gte_numb) + '_sl=' + string(sectorLength) + '_ss=' + string(sectorShift) + ..
                      '_ml=' + string(modelLength) + '_fi=' + string(forecastInterval);

steadyFileName = str_currCalcIdentif + '.' + ext_steady; // the steady mode points full file name
resTxtFileName = str_currCalcIdentif + '.' + ext_out_txt; // the text results file name full file name
//=============================================================================================================================

// CALCULATIONS
if importSteady
  [reg_all, dtm_all] = importSteadyPoints(path_steady, sep, steadyFileName);
else
  [reg_all, dtm_all] = calcSteadyPoints();
end

count_steadyModes = length(reg_all); // quantity of the all obtained steady modes points
count_initCharsPnts = length(Ngte_init); // quantity of points in initial characteristics (values in every array)

//  Define the power (N) values
if diag_sys == GTE_OIL
  N_all = p2ToPower(reg_all, count_steadyModes, p2_init, Ngte_init);
else
  N_all = ngvToPower(reg_all, Nnom, ngv_nom);
end

//  Sort Ngte parameters array values in growing order and corresponding interchange of placements values of the "dtm" parameter arrays
[N_all, dtm_all_sort] = sortByX(N_all, dtm_all);

//  Steady mode points approximation for obtaining the results characteristics
dtm_apr = approximation(N_all, dtm_all_sort, Ngte_init, polynPow, count_initCharsPnts, count_dtmParams);

//-----------------------------------
// Calc variance of the normalized steady mode points
dtm_all_apr = [];
for j = 1 : count_steadyModes
  for i = 1 : count_dtmParams
    dtm_all_apr(j, i) = interExtraPolation(Ngte_init, dtm_apr(:, i), N_all(j));
  end
end

dtm_all_dev = dtm_all_sort - dtm_all_apr; // deviation steady mode points from the approximation line
printf("variance = %f\n", variance(dtm_all_dev));

//-------------------------------------
printf("[INFO]: Characteristics was defined. Steady mode points quantity: %i\n", count_steadyModes);
//=============================================================================================================================

//      SHOW RESULTS
//  Show results processing
str_datetime = getDateTimeString();

//  Graphics plot
if plotGraphs
  plotResults(N_all, dtm_all_sort, Ngte_init, dtm_apr, ..
              count_dtmParams, plotGraphsSameWin, dt_names, 'Ngte', str_datetime);
end

//    SAVE RESULTS
// Export steady mode points in a text file
if exportSteadyPoints & ~importSteady
  saveSteadyPoints(path_steady, sep, steadyFileName, reg_all, dtm_all);
end

// Export results plots in graphics files
if exportResToImgFiles
  saveResToGraphicFiles(path_resImage, sep, str_currCalcIdentif, ext_out_images);
end

//  Save results in a text file
if exportResToTxtFile
  saveResToTextFile(path_resTxt, sep, resTxtFileName, str_currCalcIdentif, str_datetime, ..
                    polynPow, count_dtmParams, count_initCharsPnts, ..
                    Ngte_init, 'Ngte', dtm_apr, dt_names);
end

//  Show evaluating time in a console
dT = toc();
printf("\n[INFO]: Evaluating time: %i min  %4.1f sec\n", int(dT / 60), dT - int(dT / 60) * 60);

printf("**********************\n");
printf("* FINISH application *\n");
printf("**********************\n");

