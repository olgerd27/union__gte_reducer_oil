//  PERFORMING A DEFINITION OF GTE's OIL's SYSTEMS CHARACTERISTICS

tic;
xdel(winsid());
clear;
printf("*********************\n");
printf("* START application *\n");
printf("*********************\n");

//  INITIAL DATA
sheep = 2;                 // Sheep number
board = 2;                 // Board: 1 - right, 2 - left
variant = 1;               // Variant of the characteristics (for binding results in txt and image types)

sectorLength = 900;        // Length of splitting sectors
sectorShift = 300;         // Shift of sector sectorLength in every main cycle iteration
UGt_strange = 10;          // Settings Gt strange for defining the steady modes of GTE's work
Un2_xx = 5500;             // Settings XX by n2 parameter for defining the steady modes of GTE's work

exportCharactersToTxt = 0; // Export the oil's characteristics points values in a text file: 1 - perform, 0 - don't perform
exportCharactersToPNG = 0; // Export plotted the oil's characters points values in a graphics files with "png" extension:
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

path_archives = 'D:\work\GTA_M56\Archivs\sheep_2\bort_2-left_all';
ext_archive = 'txt';

// Names of oil's temperatures
t_name(index_gte_in) = 'tm_gte_in';
t_name(index_per) = 'tm_per';
t_name(index_tkvd) = 'tm_tkvd';
t_name(index_tnd) = 'tm_tnd';
t_name(index_tv) = 'tm_tv';
t_name(index_gte_out) = 'tm_gte_out';

// Paths for results saving
path_res = "D:\work\GTA_M56\documentation_preparing\programm_methods_initial_data\apps\gte_oil\out"; // common result path
path_resTxt = "txt"; // path for saving the text result
path_resImage = "images"; // path for saving the images
//=============================================================================================================================

// INITIALIZATION
paramsCount = length(params_indexes);
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
  archiveData = fscanfMat(path_archives + '\' + filesArchive(fileIndex) + '.' + ext_archive);
  
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
  for i = 1 : paramsCount
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
    for t = 1 : paramsCount
      tm_curr(:, t) = tm(from : to, t);
    end
    Ixx_curr = Ixx(from : to);
    
    Ixx_currP(j) = j * sectorShift;

    Gt_strange(j) = strange(Gt_curr1);

    // average values
    p2_avrg(j) = median(p2_curr);
    n2_avrg(j) = median(n2_curr);
    Gt_avrg(j) = median(Gt_curr2);
    for t = 1 : paramsCount
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
      for t = 1 : paramsCount
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
  for t = 1 : paramsCount - 1
    dtm_steady(:, t) = tm_steady(:, t + 1) - tm_steady(:, 1);
  end

  // Check existence the "bad" points, that is far from others points
  if length(find(dtm_steady > 80 | dtm_steady < 0))
    scf(20); xgrid;
    plot2d(p2);
    plot2d(arrayNumber, p2_steady, -6); title('Archive #' + string(fileIndex) + ' = ' + filesArchive(fileIndex), 'fontsize', 4);
    scf(21); xgrid;
    plot2d(tm);
    for i = 1 : paramsCount
      plot2d(arrayNumber, tm_steady(:, i), -6);  
    end
  end

  // Save the values of the steady mode points over all archives for further processing out of the main cycle
  p2_steadyAll = [p2_steadyAll; p2_steady];
  dtm_steadyAll = [dtm_steadyAll; dtm_steady];
  
  printf("[INFO]: Archive #%i: %s, points quantity = %i\n", fileIndex, filesArchive(fileIndex), steadyIndex);
end

countSteadyModes = length(p2_steadyAll); // quantity of the all obtained steady modes points
countInitChars = length(Ngte_init); // quantity of values in a initial characteristics

//  Define the Ngte values, using initial characteristics Ngte = f(P22)
for i = 1 : countSteadyModes
  Ngte_all(i) = sp5(p22_init, Ngte_init, p2_steadyAll(i));
end

//  Sort Ngte parameters array values in growing order and corresponding interchange of placements values of the "dtm" parameter arrays
[Ngte_steadyAllSort, dtm_steadyAllSort] = sortX(Ngte_all, dtm_steadyAll);

//  Steady mode points approximation for obtaining the oil's characteristics
//  Initialization of arrays, that will be use for inserting approximate values of characteristics line
//dtm_apr_init = approximation(Ngte_init, Ngte_steadyAllSort, polynPow);

dtm_apr_init(countInitChars, paramsCount - 1) = 0;
coef(max(polynPow) + 1, paramsCount - 1) = 0;
for j = 1 : paramsCount - 1
  coef(1 : polynPow(j) + 1, j) = approkn(Ngte_steadyAllSort, dtm_steadyAllSort(:, j), polynPow(j)); // define the values of polynomal coefficients
  for i = 1 : polynPow(j) + 1
    dtm_apr_init(:, j) = dtm_apr_init(:, j) + coef(i, j) * Ngte_init ^ (polynPow(j) + 1 - i);  // getting the results oil's characteristics
  end
end

printf("\n");
printf("[INFO]: Characteristics was defined. Steady mode points quantity: %i\n", countSteadyModes);
//=============================================================================================================================

//      SHOW RESULTS
//  Show results in console
//printf("RESULTS:\n");
strTitle = "sheep = " + string(sheep) + ", board = " + string(board) + "\nsectorLength = " + string(sectorLength) + ..
           ", sectorShift = " + string(sectorShift) + "\nvariant = " + string(variant) + "\n";
//printf("\n%s\n", strTitle);

//  Graphics plot
for t = 1 : paramsCount - 1
  strImagesTitle = 'd' + t_name(t) + ' = f(Ngte)';
  hWin = scf(t); plot2d(Ngte_steadyAllSort, dtm_steadyAllSort(:, t), -9); xgrid; 
  title(strImagesTitle + ',  variant = ' + string(variant), 'fontsize', 4); // steady mode points
  hWin.figure_name = strImagesTitle;

  //  Plot the approximate line
  plot2d(Ngte_init, dtm_apr_init(:, t), 15); e = gce(); e.children.thickness = 2;
  // Plot the approximate lines points
  plot2d(Ngte_init, dtm_apr_init(:, t), -14);
end

//    SAVE RESULTS
// Export results plots in graphics files
if exportCharactersToPNG == 1
  // dirs create/delete operations
  // go into the common directory for saving all images
  path_resImage = path_res + '\' + path_resImage;
  if(~isdir(path_resImage) & ~createdir(path_resImage))
    printf("[ERROR]: Cannot create the common dir for saving result images:\n\t%s\n", path_resImage);
  end
  
  // go into the special directory for current calculation
  path_resImage = path_resImage + '\' + 'sheep' + string(sheep) + '_board' + string(board);
  if isdir(path_resImage)
    removedir(path_resImage);
  end
  if(~createdir(path_resImage))
    printf("[ERROR]: Cannot create the dir for saving result images: %s\n", path_resImage);
  end
  
  // save image data
  printf("[INFO]: Export the results plotted images of the GTE''s oil''s characteristics to the files:\n");
  figureIDs = winsid();
  for i = 1 : length(figureIDs)
    h = scf(figureIDs(i));
    exportFileName = h.figure_name + '.png';
    exportFileName = strsubst(exportFileName, ' ', ''); // spaces deleting
    filePathName = path_resImage + '\' + exportFileName;
    xs2png(h, filePathName);
    printf("\t%s\n", filePathName);
  end
end

//  Save results in a text file
if exportCharactersToTxt == 1
  // initial data
  path_resTxt = path_res + '\' + path_resTxt;
  file_resTxt = 'sheep' + string(sheep) + '_board' + string(board) + '.rez';
  pathFile_resTxt = path_resTxt + '\' + file_resTxt;
  
  if ~isdir(path_resTxt)
    createdir(path_resTxt);
  else
    mdelete(pathFile_resTxt); // delete file with the same name as a current file
  end
  
  f = mopen(pathFile_resTxt, 'w');
  mfprintf(f, strTitle + "powers = ");
  for i = 1 : paramsCount - 1
    mfprintf(f, "%i  ", polynPow(i));
  end
  
  // Write the parameters names
  mfprintf(f, "\n\t\tNgte  ");
  for i = 1 : paramsCount - 1
    mfprintf(f, "d%s\t\t", t_name(i));
  end
  mfprintf(f, "\n");
  
  // Write values
  for i = 1 : countInitChars
    mfprintf(f, "%8.2f\t\t", Ngte_init(i));
    for j = 1 : paramsCount - 1
      mfprintf(f, "%5.2f\t\t\t", dtm_apr_init(i, j));
    end
    mfprintf(f, "\n");
  end
  mclose(f);
  
  // show info message
  printf("[INFO]: The results points values of the GTE''s oil''s characteristics was exported to the text file:\n");
  printf("\t%s\n", pathFile_resTxt);
end

//  Show evaluating time in a console
printf("**********************\n");
printf("* FINISH application *\n");
printf("**********************\n");

dT = toc();
printf("\n[INFO]: Evaluating time: %i min  %4.1f sec\n", int(dT / 60), dT - int(dT / 60) * 60);
printf("========================================================================================================\n\n");

