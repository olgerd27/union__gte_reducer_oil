//  PERFORMING A DEFINITION OF GTE's and OIL's SYSTEMS CHARACTERISTICS

tic;
xdel(winsid());
clear;
stacksize(5e7);
printf("*********************\n");
printf("* START application *\n");
printf("*********************\n");

//  INITIAL DATA
diag_sys = 2;             // PKSTD diagnostics system: 1 - GTE oil's system, 2 - reducer oil's system
sheep = 2;                // Sheep number
board = 2;                // Board: 1 - right, 2 - left

sectorLength = 900;       // Length of splitting sectors
sectorShift = 300;        // Shift of sector with length sectorLength in every main cycle iteration
if diag_sys == 1
  UGt_strange = 10;       // Settings Gt strange for defining the steady modes of GTE's work
  Un2_xx = 5500;          // Settings XX by n2 parameter for defining the steady modes of GTE's work
else
  Ungv_strange = 5;       // Settings ngv strange for defining the steady modes of reducer's work
  Ungv_min = 40;          // Settings ngv min value for defining the steady modes of reducer's work

  Nnom = 20020;           // the reducer power on the nominal reducer's work mode
  ngv_nom = 240;          // rotation speed of the reducer outlet shaft on the nominal reducer's work mode
end

Udtm_valid_min = 0;       // Settings dtm min value for defining the invalid points
Udtm_valid_max = 80;      // Settings dtm max value for defining the invalid points


exportResToTxtFile  = 0;  // Export the oil's characteristics points values in a text file: 1 - perform, 0 - don't perform
exportResToImgFiles = 0;  // Export plotted the oil's characters points values in a graphics files with "png" extension:
                          // 1 - perform, 0 - don't perform

// The initial characteristics Ngte = f(p2) from a set of thermodynamics characteristics
p2_init = [1.795; 2.251; 2.702; 3.380; 4.147; 5.111; 6.322; 7.049; 7.775; 8.408;
            9.152; 9.836; 10.601; 12.189; 13.692; 15.132; 16.397; 17.560; 18.671; 19.782];
Ngte_init = [110.4; 273.4; 434.4; 676.5; 950.2; 1294.5; 2052.5; 2643.4; 3234.2; 3749.1;
             4353.8; 5098.5; 5958.6; 7823; 9869.7; 11981; 14013.4; 16015.3; 18017.2; 20019.1];

// Indexes of the oil's parameters
index_in = 0; index_out = 0; // initialization
if diag_sys == 1
  // TODO: Maybe delete all indexes and params_indexes here and in "else" section too... it is not used anywhere but INITIAL DATA section
  index_in = 1;  index_per = 2;  index_tkvd = 3;  index_tnd = 4;  index_tv = 5;  index_out = 6;
  params_indexes = [index_in; index_per; index_tkvd; index_tnd; index_tv; index_out];
else
  index_in = 1;      index_out = 2;
  index_tz01a =  3;  index_tz01b =  4;  index_tz02a =  5;  index_tz02b =  6;  index_tz02c =  7;  index_tz02d =  8;
  index_tz03a =  9;  index_tz03b = 10;  index_tz03c = 11;  index_tz03d = 12;  index_tz04a = 13;  index_tz04b = 14;
  index_tz05a = 15;  index_tz05b = 16;  index_tz06a = 17;  index_tz06b = 18;  index_tz07a = 19;  index_tz07b = 20;
  index_tz08a = 21;  index_tz08b = 22;  index_tz09a = 23;  index_tz09b = 24;  index_tz10a = 25;  index_tz10b = 26;
  index_tz11a = 27;  index_tz11b = 28;
  params_indexes = [index_in; index_out; 
    index_tz01a; index_tz01b; index_tz02a; index_tz02b; index_tz02c; index_tz02d; index_tz03a; index_tz03b; 
    index_tz03c; index_tz03d; index_tz04a; index_tz04b; index_tz05a; index_tz05b; index_tz06a; index_tz06b;
    index_tz07a; index_tz07b; index_tz08a; index_tz08b; index_tz09a; index_tz09b; index_tz10a; index_tz10b;
    index_tz11a; index_tz11b];
end

// Polynomial powers that describe the oil's characteristics dtX = f(Ngte)
// TODO: maybe use simply digit indexes instead of variables. This make the code more simply
if diag_sys == 1
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

// Setting archives path, names and extension
filesArchive = [// 1 etap PI
                'mo_2013_1_3_23_0_0';   'mo_2013_1_3_9_52_48';  'mo_2013_1_4_9_34_5'; 'mo_2013_1_5_9_27_47';
                'mo_2013_1_8_11_10_49'; 'mo_2013_1_8_16_25_21';
                // 2 etap PI + PSI
                'mo_2013_2_14_10_1_0';  'mo_2013_2_25_13_24_20';  'mo_2013_2_25_18_44_48'; 'mo_2013_2_25_9_58_3';
                'mo_2013_2_26_9_28_29'; 'mo_2013_3_11_12_48_00'];

path_archives = "/media/oleg/users/Oleg/work_zm/export/GTA_M56/Archivs/sheep2/sheep2_bort2/bort_2-all";
//path_archives = "D:\work\GTA_M56\Archivs\sheep_2\bort_2-left_all";

// Names of oil's temperatures
if diag_sys == 1
  t_name(index_in) = 'tm_gte_in';  t_name(index_per) = 'tm_per';
  t_name(index_tkvd) = 'tm_tkvd';  t_name(index_tnd) = 'tm_tnd';
  t_name(index_tv) = 'tm_tv';      t_name(index_out) = 'tm_gte_out';
else
  t_name(index_in) = 'tm_red_in';  t_name(index_out) = 'tm_red_out';
  t_name(index_tz01a) = 'tz01a';   t_name(index_tz01b) = 'tz01b';
  t_name(index_tz02a) = 'tz02a';   t_name(index_tz02b) = 'tz02b';
  t_name(index_tz02c) = 'tz02c';   t_name(index_tz02d) = 'tz02d';
  t_name(index_tz03a) = 'tz03a';   t_name(index_tz03b) = 'tz03b';
  t_name(index_tz03c) = 'tz03c';   t_name(index_tz03d) = 'tz03d';
  t_name(index_tz04a) = 'tz04a';   t_name(index_tz04b) = 'tz04b';
  t_name(index_tz05a) = 'tz05a';   t_name(index_tz05b) = 'tz05b';
  t_name(index_tz06a) = 'tz06a';   t_name(index_tz06b) = 'tz06b';
  t_name(index_tz07a) = 'tz07a';   t_name(index_tz07b) = 'tz07b';
  t_name(index_tz08a) = 'tz08a';   t_name(index_tz08b) = 'tz08b';
  t_name(index_tz09a) = 'tz09a';   t_name(index_tz09b) = 'tz09b';
  t_name(index_tz10a) = 'tz10a';   t_name(index_tz10b) = 'tz10b';
  t_name(index_tz11a) = 'tz11a';   t_name(index_tz11b) = 'tz11b';
end

// Paths for results saving
path_res = "~/Programming/scilab/projects/union__gte_reducer_oil/out";
//path_res = "D:\work\GTA_M56\documentation_preparing\programm_methods_initial_data\apps\union__gte_reducer_oil\out";
//=============================================================================================================================

// LOADING additional files with functions
path_sourceFiles = "~/Programming/scilab/projects/union__gte_reducer_oil/src";
//path_sourceFiles = "D:\work\GTA_M56\documentation_preparing\programm_methods_initial_data\apps\union__gte_reducer_oil\src";
names_sourceFiles = [
                    "add_functions.sci";
                    "in_out_functions.sci";
                    "special_functions.sci";
                    ];
for i = 1 : size(names_sourceFiles, 'r')
  exec(path_sourceFiles + '/' + names_sourceFiles(i)); // loading functionality from the files
end
//=============================================================================================================================

// INITIALIZATION
if (index_in == 0) | (index_out == 0)
  printf("[INFO]: The ""index_in"" or ""index_out"" is not defined in the INITIAL DATA section\n");
  return;
end
count_tmParams = length(params_indexes); // quantity of the temperature
count_dtmParams = count_tmParams - 1; // quantity of the temperatures delta's
for i = 1 : count_dtmParams
  dt_name(i) = 'd' + t_name(i + 1); // names of the temperatures delta's
end
ext_archive = 'txt'; // archives files extention
// Arrays for storing steady mode points values
// TODO: is this need??
reg_all = [];
dtm_all = [];
// Structure for storing parameters data, readed from an archive file
if diag_sys == 1
  index_t0 = 1;  index_Gt = 2;  index_reg = 3;  index_n2 = 4;  index_tm = 5; // parameters indexes
  params(index_t0) = struct('name', 't0', 'archIndexStart', 16, 'archIndexEnd', 18, 'data', []);
  params(index_Gt) = struct('name', 'Gt', 'archIndexStart', 120, 'archIndexEnd', 120, 'data', []);
  params(index_reg) = struct('name', 'p2', 'archIndexStart', 38, 'archIndexEnd', 38, 'data', []);
  params(index_n2) = struct('name', 'n2', 'archIndexStart', 10, 'archIndexEnd', 10, 'data', []);
  params(index_tm) = struct('name', t_name, 'archIndexStart', 47, 'archIndexEnd', 52, 'data', []);
else
  index_reg = 1;  index_tm = 2; // parameters indexes
  params(index_reg) = struct('name', 'ngv', 'archIndexStart', 14, 'archIndexEnd', 14, 'data', []);
  params(index_tm) = struct('name', t_name, 'archIndexStart', 53, 'archIndexEnd', 80, 'data', []);
end
colors = [1, 2, 3, 5, 19, 16, 27, 22, 13, 6, 9, 32, 28, 21, 25, 23, 26, 17];
//=============================================================================================================================

// MAIN CYCLE
for fileIndex = 1 : 3//size(filesArchive, 'r') // TODO: correct cycle quantity after refactoring
  // Deleting results that was obtained in the previous main cycle iteration
  clear reg_steady; clear tm_steady; clear dtm_steady; clear arrayNumber_steady;
  
  //  READING an archive file
  params = readParametersData(path_archives + '/', filesArchive(fileIndex), '.' + ext_archive, params);

  // Getting the parameters arrays with reduction to a normal atmospheric condition
  // This parameters is take out of the if-else, because name of its index variable for the both diagnostic systems is equal
  reg = params(index_reg).data; // the regime parameter values
  tm = params(index_tm).data;
  if diag_sys == 1
    t0 = mean(params(index_t0).data, 'c');
    alpha = sqrt(288 ./ (t0 + 273)); // alpha coefficient for reductions parameters to normal atmospheric conditions
    Gt = params(index_Gt).data;
    reg = reg * 10.2;  // there are p2 parameter values with conversion its from MPa to kg/cm2
    n2 = params(index_n2).data .* alpha; // reduction to normal conditions
  end

  //  Splitting parameters arrays to sectors with defining the average values and strange
  steadyIndex = 0;
  kk = ceil(sectorLength / sectorShift); // the coefficient for correct definition the iterations quantity and correct shifting buffer with length "sectorShift" along the full length of archive
  arrSize = length(reg);
  for j = kk : int(arrSize / sectorShift)
    from = sectorShift * (j - 1) + 1;
    to = sectorShift * j;
    arrayNumber = j * sectorShift;

    // split arrays, calc strange or average values of parameters and define steady mode
    isSteadyMode = %F;
    if diag_sys == 1
      Gt_strange = strange(Gt(to - sectorLength + 1 : to)); // there are splitting and strange value calculation
      n2_avrg = median(n2(from : to));
      Gt_avrg = median(Gt(from : to));
      isSteadyMode = (Gt_strange <= UGt_strange) & (n2_avrg > Un2_xx) & (Gt_avrg > 0); // define steady mode
    else
      reg_strange = strange(reg(to - sectorLength + 1 : to));
      reg_avrg = median(reg(from : to));
      isSteadyMode = (reg_strange <= Ungv_strange) & (reg_avrg > Ungv_min); // define steady mode
    end

    // Processing the steady modes of GTE's work points values
    if isSteadyMode
      steadyIndex = steadyIndex + 1;
      if diag_sys == 1
        reg_steady(steadyIndex) = median(reg(from : to));
      else
        reg_steady(steadyIndex) = reg_avrg; // already calculated
      end
      for t = 1 : count_tmParams
        tm_steady(steadyIndex, t) = median(tm(from : to, t));
      end
      arrayNumber_steady(steadyIndex) = arrayNumber;
    end
  end
  
  //  Check if in the current archive doesn't exist the steady mode points
  if steadyIndex == 0
    printf("[ERROR]: Steady mode points not found: archive #%i = %s, sectorLength = %i\n", fileIndex, filesArchive(fileIndex), sectorLength);
    printf("Continue? (1 - yes, 2 - no)\n");
    key = scanf("%i");
    if key == 1
      continue;
    else
      scf(1); xgrid; title('Steady mode points not found. ' + params(index_reg).name + ' = f(time)', 'fontsize', 4);
      plot2d(reg);  e = gce(); e.children.thickness = 2;
      printf("--------------------------------------------------------------------------------------------------------\n\n");
      return;
    end
  end

  //  Define temperature drop of oil
  for t = 1 : count_dtmParams
    dtm_steady(:, t) = tm_steady(:, t + 1) - tm_steady(:, index_in);
  end

  // TODO: check this code, maybe there are need to refactor this
  // Check existence the "bad", invalid points, that is far from others points
  ind_invalidValues = find(dtm_steady > Udtm_valid_max | dtm_steady < Udtm_valid_min);
  count_invalidPoints = length(ind_invalidValues);
  if count_invalidPoints
    [cols_invalid_dt, rows_invalid] = calcInvalidValuePos(ind_invalidValues, size(reg_steady, 'r'));
    rows_invalid_u = unique(rows_invalid);
    str_archiveNumberName = 'archive #' + string(fileIndex) + ': ' + filesArchive(fileIndex)';

    printf("[ERROR]: There was found %i invalid steady mode point(s) in the %s\n", count_invalidPoints, str_archiveNumberName);
    for i = 1 : count_invalidPoints
      printf("\tpoint #%i: parameter = ''%s'', number = %i\n", i, dt_name(cols_invalid_dt(i)), rows_invalid(i));
    end
    
    printf("Continue? (1 - yes, 2 - no)\n");
    key = scanf("%i");
    if key == 1
      // delete rows with invalid point(-s) for getting arrays with only valid points
      reg_steady(rows_invalid_u, :) = [];
      dtm_steady(rows_invalid_u, :) = [];
      continue;
    else
      plotInvalidArchive( reg, tm, reg_steady, tm_steady, arrayNumber_steady, ..
                                       cols_invalid_dt, rows_invalid, rows_invalid_u, ..
                                       index_in, params(index_reg).name, t_name, ..
                                       str_archiveNumberName, colors );
      printf("--------------------------------------------------------------------------------------------------------\n\n");
      return;
    end
  end

  // Save the values of the steady mode points over all archives for further processing out of the main cycle
  reg_all = [reg_all; reg_steady];
  dtm_all = [dtm_all; dtm_steady];
    
  printf("[INFO]: Archive #%i: ""%s"", points quantity = %i\n", fileIndex, filesArchive(fileIndex), steadyIndex);
end

count_steadyModes = length(reg_all); // quantity of the all obtained steady modes points
count_initCharsPnts = length(Ngte_init); // quantity of points in initial characteristics (values in every array)

//  Define the power (N) values
if diag_sys == 1
  N_all = p2ToPower(reg_all, count_steadyModes, p2_init, Ngte_init);
else
  N_all = ngvToPower(reg_all, Nnom, ngv_nom);
end

//  Sort Ngte parameters array values in growing order and corresponding interchange of placements values of the "dtm" parameter arrays
[N_all, dtm_all] = sortByX(N_all, dtm_all);

//  Steady mode points approximation for obtaining the results characteristics
dtm_apr = approximation(N_all, dtm_all, Ngte_init, polynPow, count_initCharsPnts, count_dtmParams);

printf("[INFO]: Characteristics was defined. Steady mode points quantity: %i\n", count_steadyModes);
//=============================================================================================================================

//      SHOW RESULTS
// CONSTANT INITIAL DATA
// Relative path's for saving results. 
// Structure of the results tree directories is constant and therefore this paths variables moved from a Initial data block at here
path_resGTERltv = "gte"; // relative path for saving the GTE's results characteristics
path_resReducerRltv = "reducer"; // relative path for saving the reducer's results characteristics
path_resTxtRltv = "txt"; // relative path for saving the text result
path_resImageRltv = "images"; // relative path for saving the images
// Results files exensions
ext_images = 'png';
ext_txt = 'rez';
// Identification current calculation - sheep+board
currentCalcIdentif = 'sheep' + string(sheep) + '_board' + string(board);

// Define the result relative path in accordance with type of current PKSTD diagnostics system
if diag_sys == 1
  path_resDiagSysRltv = path_resGTERltv;
else
  path_resDiagSysRltv = path_resReducerRltv;
end

//  Show results processing
str_datetime = getDateTimeString();
strTitle = currentCalcIdentif + '\nsectorLength = ' + string(sectorLength) + ..
           ', sectorShift = ' + string(sectorShift) + '\ndate_time = ' + str_datetime + '\n';
//printf("\n%s\n", strTitle);

//  Graphics plot
windowNumber = max(winsid()) + 1;
for t = 1 : count_dtmParams
  strImagesTitle = dt_name(t) + ' = f(Ngte)';
  hWin = scf(windowNumber + t); plot2d(N_all, dtm_all(:, t), -9); xgrid; 
  title(strImagesTitle + ',  ' + str_datetime, 'fontsize', 4); // steady mode points
  hWin.figure_name = strImagesTitle;

  //  Plot the approximate line
  plot2d(Ngte_init, dtm_apr(:, t), 15); e = gce(); e.children.thickness = 2;
  // Plot the approximate lines points
  plot2d(Ngte_init, dtm_apr(:, t), -14);
end

//    SAVE RESULTS
// Export results plots in graphics files
if exportResToImgFiles == 1
  saveResToGraphicFiles(path_res + '/' + path_resDiagSysRltv + '/' + path_resImageRltv, ..
                        currentCalcIdentif, ext_images);
end

//  Save results in a text file
if exportResToTxtFile == 1
  saveResToTextFile(path_res + '/' + path_resDiagSysRltv + '/' + path_resTxtRltv, ..
                    currentCalcIdentif + '.' + ext_txt, strTitle, polynPow, ..
                    count_dtmParams, count_initCharsPnts, ..
                    Ngte_init, 'Ngte', dtm_apr, dt_name);
end

//  Show evaluating time in a console
dT = toc();
printf("\n[INFO]: Evaluating time: %i min  %4.1f sec\n", int(dT / 60), dT - int(dT / 60) * 60);

printf("**********************\n");
printf("* FINISH application *\n");
printf("**********************\n");

