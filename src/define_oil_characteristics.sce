//  PERFORMING A DEFINITION OF GTE's and OIL's SYSTEMS CHARACTERISTICS

tic;
xdel(winsid());
clear;
printf("*********************\n");
printf("* START application *\n");
printf("*********************\n");

// LOADING additional files with functions
//path_sourceFiles = "~/Programming/scilab/projects/union__gte_reducer_oil/src";
path_sourceFiles = "D:\work\GTA_M56\documentation_preparing\programm_methods_initial_data\apps\union__gte_reducer_oil\src";
names_sourceFiles = [
                    "add_functions.sci";
                    ];
for i = 1 : size(names_sourceFiles, 'r')
  exec(path_sourceFiles + '/' + names_sourceFiles(i)); // loading functionality from the files
end

//  INITIAL DATA
diag_sys = 1;             // PKSTD system type: 1 - GTE oil's system, 2 - reducer oil's system
sheep = 2;                // Sheep number
board = 2;                // Board: 1 - right, 2 - left

sectorLength = 900;       // Length of splitting sectors
sectorShift = 300;        // Shift of sector sectorLength in every main cycle iteration
UGt_strange = 10;         // Settings Gt strange for defining the steady modes of GTE's work
Un2_xx = 5500;            // Settings XX by n2 parameter for defining the steady modes of GTE's work

exportResToTxtFile = 0;   // Export the oil's characteristics points values in a text file: 1 - perform, 0 - don't perform
exportResToImgFiles = 1;  // Export plotted the oil's characters points values in a graphics files with "png" extension:
                          // 1 - perform, 0 - don't perform

// The initial characteristics Ngte = f(P22) from a set of thermodynamics characteristics
p22_init = [1.795; 2.251; 2.702; 3.380; 4.147; 5.111; 6.322; 7.049; 7.775; 8.408;
            9.152; 9.836; 10.601; 12.189; 13.692; 15.132; 16.397; 17.560; 18.671; 19.782];
Ngte_init = [110.4; 273.4; 434.4; 676.5; 950.2; 1294.5; 2052.5; 2643.4; 3234.2; 3749.1;
             4353.8; 5098.5; 5958.6; 7823; 9869.7; 11981; 14013.4; 16015.3; 18017.2; 20019.1];

// Indexes of the oil's parameters
index_gte_in = 1; 
index_per = 2;  index_tkvd = 3;  index_tnd = 4;  index_tv = 5;  index_gte_out = 6;
params_indexes = [index_gte_in; index_per; index_tkvd; index_tnd; index_tv; index_gte_out];

// Polynomial powers that describe the oil's characteristics dtm_X = f(Ngte)
polynPow(index_per - index_gte_in) = 2;     // dtm_per
polynPow(index_tkvd - index_gte_in) = 2;    // dtm_tkvd
polynPow(index_tnd - index_gte_in) = 2;     // dtm_tnd
polynPow(index_tv - index_gte_in) = 2;      // dtm_tv
polynPow(index_gte_out - index_gte_in) = 2; // dtm_gte_out

// Setting archives path, names and extension
filesArchive = [// 1 etap PI
                'mo_2013_1_3_23_0_0';   'mo_2013_1_3_9_52_48';  'mo_2013_1_4_9_34_5'; 'mo_2013_1_5_9_27_47';
                'mo_2013_1_8_11_10_49'; 'mo_2013_1_8_16_25_21';
                // 2 etap PI + PSI
                'mo_2013_2_14_10_1_0';  'mo_2013_2_25_13_24_20';  'mo_2013_2_25_18_44_48'; 'mo_2013_2_25_9_58_3';
                'mo_2013_2_26_9_28_29'; 'mo_2013_3_11_12_48_00'];

//path_archives = '/media/users/Oleg/work_zorya-mashproekt/export20130702/GTA_M56/Archivs/sheep2/sheep2_bort2/bort_2-all';
path_archives = "D:\work\GTA_M56\Archivs\sheep_2\bort_2-left_all";
ext_archive = 'txt';

// Names of oil's temperatures
t_name(index_gte_in) = 'tm_gte_in';
t_name(index_per) = 'tm_per';
t_name(index_tkvd) = 'tm_tkvd';
t_name(index_tnd) = 'tm_tnd';
t_name(index_tv) = 'tm_tv';
t_name(index_gte_out) = 'tm_gte_out';

// Paths for results saving
//path_res = "/media/users/Oleg/work_zorya-mashproekt/export20130702/GTA_M56/documentation_preparing/programm_methods_initial_data/apps/union__gte_reducer_oil/out";
path_res = "D:\work\GTA_M56\documentation_preparing\programm_methods_initial_data\apps\union__gte_reducer_oil\out";
//=============================================================================================================================

// INITIALIZATION
count_tmParams = length(params_indexes);
count_dtmParams = count_tmParams - 1;
// Arrays for storing steady mode points values
p2_steadyAll = [];
dtm_steadyAll = [];

// MAIN CYCLE
for fileIndex = 1 : size(filesArchive, 'r')
  // Deleting results that was obtained in the previous main cycle iteration
  clear archiveData; clear t0; clear alfa; clear Gt; clear n2; clear tm;
  clear Ixx_currP; clear Gt_strange;
  clear p2_avrg; clear n2_avrg; clear Gt_avrg; clear tm_avrg;
  clear p2_steady; clear n2_steady; clear tm_steady; clear arrayNumber;
  clear dtm_steady;
  
  //  READING an archive file
  archiveData = fscanfMat(path_archives + '/' + filesArchive(fileIndex) + '.' + ext_archive);
  
  // Reading data from a bad archives
  if filesArchive(fileIndex) == 'mo_2012_11_24_8_55_21'
    // In this archive was breaked one sensor. There are need to use only part of archive, when values of this sensors was valid
    archiveData = archiveData(1 : 8750, :);
  end

  // Getting the parameters arrays with reduction to a normal atmospheric condition
  t0 = mean(archiveData(:, 16:18), 'c');
  alfa = sqrt(288 ./ (t0 + 273));
  Gt = archiveData(:, 120);
  p2 = archiveData(:, 38) * 10.2;  // for getting a regime parameters values (Ngte), there are needing to read P22 values. Readed values P22 converting from MPa to kg/cm2
  n2 = archiveData(:, 10) .* alfa;
  for i = 1 : count_tmParams
    tm(:, i) = archiveData(:, 46 + i);
  end

  IN = length(Gt);
  Ixx = [1 : IN];
  
  //  Splitting of parameters arrays to sectors with defining the average values and strange
  kk = ceil(sectorLength / sectorShift); // the coefficient for correct definition the iterations quantity and correct moving buffer with length "sectorShift" secs by all archive
  for j = kk : int(IN / sectorShift)
    from = sectorShift * (j - 1) + 1;
    to = sectorShift * j;
    
    // splitting arrays
    Gt_curr1 = Gt(to - sectorLength + 1 : to);
    Gt_curr2 = Gt(from : to);
    p2_curr = p2(from : to);
    n2_curr = n2(from : to);
    for t = 1 : count_tmParams
      tm_curr(:, t) = tm(from : to, t);
    end
    Ixx_curr = Ixx(from : to);
    
    Ixx_currP(j) = j * sectorShift;

    Gt_strange(j) = strange(Gt_curr1);

    // average values
    p2_avrg(j) = median(p2_curr);
    n2_avrg(j) = median(n2_curr);
    Gt_avrg(j) = median(Gt_curr2);
    for t = 1 : count_tmParams
      tm_avrg(j, t) = median(tm_curr(:, t));
    end
  end
     
  // Define the steady modes of GTE's work on sectors with length "sectorLength" seconds
  steadyIndex = 0;
  for i = 1 : int(IN / sectorShift)
    if (Gt_strange(i) <= 10) & (n2_avrg(i) > Un2_xx) & (Gt_avrg(i) > 0)
      steadyIndex = steadyIndex + 1;
      p2_steady(steadyIndex) = p2_avrg(i);
      n2_steady(steadyIndex) = n2_avrg(i);
      for t = 1 : count_tmParams
        tm_steady(steadyIndex, t) = tm_avrg(i, t);
      end
      arrayNumber(steadyIndex) = Ixx_currP(i);
    end
  end
  
  //  Check if in the current archive doesn't exist the steady mode points
  if steadyIndex == 0
    scf(100); xgrid; plot2d(p2); title('Steady mode points not found', 'fontsize', 4);
    printf("Archive #%i = %s, steady mode points not found, sectorLength = %i\n", fileIndex, filesArchive(fileIndex), sectorLength);
    printf("--------------------------------------------------------------------------------------------------------\n\n");
    return;
  end

  //  Define temperature drop of oil
  for t = 1 : count_dtmParams
    dtm_steady(:, t) = tm_steady(:, t + 1) - tm_steady(:, 1);
  end

  // Check existence the "bad" points, that is far from others points
  if length(find(dtm_steady > 80 | dtm_steady < 0))
    scf(20); xgrid;
    plot2d(p2);
    plot2d(arrayNumber, p2_steady, -6); title('Archive #' + string(fileIndex) + ' = ' + filesArchive(fileIndex), 'fontsize', 4);
    scf(21); xgrid;
    plot2d(tm);
    for i = 1 : count_tmParams
      plot2d(arrayNumber, tm_steady(:, i), -6);  
    end
  end

  // Save the values of the steady mode points over all archives for further processing out of the main cycle
  p2_steadyAll = [p2_steadyAll; p2_steady];
  dtm_steadyAll = [dtm_steadyAll; dtm_steady];
  
  printf("[INFO]: Archive #%i: ""%s"", points quantity = %i\n", fileIndex, filesArchive(fileIndex), steadyIndex);
end

count_steadyModes = length(p2_steadyAll); // quantity of the all obtained steady modes points
count_initCharsPnts = length(Ngte_init); // quantity of points in initial characteristics (values in every array)

//  Define the Ngte values, using initial characteristics Ngte = f(P22)
for i = 1 : count_steadyModes
  Ngte_all(i) = interExtraPolation(p22_init, Ngte_init, p2_steadyAll(i));
end

//  Sort Ngte parameters array values in growing order and corresponding interchange of placements values of the "dtm" parameter arrays
[Ngte_all, dtm_steadyAll] = sortByX(Ngte_all, dtm_steadyAll);

//  Steady mode points approximation for obtaining the results characteristics
dtm_apr = approximation(Ngte_all, dtm_steadyAll, Ngte_init, polynPow, count_initCharsPnts, count_dtmParams);

printf("\n");
printf("[INFO]: Characteristics was defined. Steady mode points quantity: %i\n", count_steadyModes);
//=============================================================================================================================

//      SHOW RESULTS
// Relative path's for saving results. 
// Structure of the results tree directories is constant and therefore this paths variables moved from a Initial data block at here
path_resGTERltv = "gte"; // relative path for saving the GTE's results characteristics
path_resReducerRltv = "reducer"; // relative path for saving the reducer's results characteristics
path_resTxtRltv = "txt"; // relative path for saving the text result
path_resImageRltv = "images"; // relative path for saving the images

// Define the result relative path in accordance with type of current PKSTD diagnostics system
if diag_sys == 1
  path_resDiagSysRltv = path_resGTERltv;
else
  path_resDiagSysRltv = path_resReducerRltv;
end

//  Show results in console
//printf("RESULTS:\n");
str_datetime = getDateTimeString();
strTitle = "sheep = " + string(sheep) + ", board = " + string(board) + "\nsectorLength = " + string(sectorLength) + ..
           ", sectorShift = " + string(sectorShift) + "\ndate_time = " + str_datetime + "\n";
//printf("\n%s\n", strTitle);

//  Graphics plot
for t = 1 : count_dtmParams
  strImagesTitle = 'd' + t_name(t) + ' = f(Ngte)';
  hWin = scf(t); plot2d(Ngte_all, dtm_steadyAll(:, t), -9); xgrid; 
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
  // Dirs create/delete operations
  // go into the common directory for saving all images
  path_resImage = path_res + '/' + path_resDiagSysRltv + '/' + path_resImageRltv;
  if(~isdir(path_resImage) & ~createdir(path_resImage))
    printf("[ERROR]: Cannot create the common dir for saving result images:\n\t%s\n", path_resImage);
  end
  // go into the special directory for current calculation
  path_resImage = path_resImage + '/' + 'sheep' + string(sheep) + '_board' + string(board);
  if isdir(path_resImage)
    removedir(path_resImage);
  end
  if(~createdir(path_resImage))
    printf("[ERROR]: Cannot create the dir for saving result images: %s\n", path_resImage);
  end
  
  // save image data
  printf("[INFO]: Export the results plotted images of the GTE''s oil''s characteristics to the files:\n");
  ext_images = 'png';
  figureIDs = winsid();
  for i = 1 : length(figureIDs)
    h = scf(figureIDs(i));
    exportFileName = h.figure_name + '.' + ext_images;
    exportFileName = strsubst(exportFileName, ' ', ''); // spaces deleting
    filePathName = path_resImage + '/' + exportFileName;
    xs2png(h, filePathName);
    printf("\t%s\n", filePathName);
  end
end

//  Save results in a text file
if exportResToTxtFile == 1
  // initial data
  path_resTxt = path_res + '/' + path_resTxtRltv;
  file_resTxt = 'sheep' + string(sheep) + '_board' + string(board) + '.rez';
  pathFile_resTxt = path_resTxt + '/' + file_resTxt;
  
  if ~isdir(path_resTxt)
    createdir(path_resTxt);
  else
    mdelete(pathFile_resTxt); // delete file with the same name as a current file
  end
  
  f = mopen(pathFile_resTxt, 'w');
  mfprintf(f, strTitle + "powers = ");
  for i = 1 : count_dtmParams
    mfprintf(f, "%i  ", polynPow(i));
  end
  
  // Write the parameters names
  mfprintf(f, "\n\t\tNgte  ");
  for i = 1 : count_dtmParams
    mfprintf(f, "d%s\t\t", t_name(i));
  end
  mfprintf(f, "\n");
  
  // Write values
  for i = 1 : count_initCharsPnts
    mfprintf(f, "%8.2f\t\t", Ngte_init(i));
    for j = 1 : count_dtmParams
      mfprintf(f, "%5.2f\t\t\t", dtm_apr(i, j));
    end
    mfprintf(f, "\n");
  end
  mclose(f);
  
  // show info message
  printf("[INFO]: The results points values of the GTE''s oil''s characteristics was exported to the text file:\n");
  printf("\t%s\n", pathFile_resTxt);
end

//  Show evaluating time in a console
dT = toc();
printf("\n[INFO]: Evaluating time: %i min  %4.1f sec\n", int(dT / 60), dT - int(dT / 60) * 60);

printf("**********************\n");
printf("* FINISH application *\n");
printf("**********************\n");

