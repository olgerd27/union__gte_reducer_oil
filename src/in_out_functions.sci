// Functions for performing the in-, out- operations

// ======= DATA READING =======
function [reg_all, dtm_all] = importSteadyPoints(path, fileName)
//*******************************************************
// 
//*******************************************************
  filePathName = path + '/' + fileName;
  
  if ~isfile(filePathName)
    printf("[ERROR]: cannot read the steady mode points data from the file:\n");
    printf("\t''%s''\n", filePathName);
    printf("\tThe reason: this file is not exist or it is not a file\n");
    abort;
  end

  data = fscanfMat(filePathName);
  reg_all = data(:, 1);
  dtm_all = data(:, 2 : size(data, 'c'));
endfunction

function resData = getValidData(fileName, initData)
//***************************************************************************
// Get a valid data from some archives with an invalid data.                *
// The problem may to be for example - breaking some sensor.                *
//                                                                          *
// The knowledge about invalid archives and diapasons of a valid            *
// data in it was moved from a main function to here, that make it simpler. *
//                                                                          *
// IN:  fileName - archive file name                                        *
//      initData - initial data, readed from a archive file                 *
// OUT: resData - result data                                               *
//***************************************************************************
  // initialization
  from = 1;
  to = size(initData, 'r');

  // checking archive data
  if diag_sys == 1
    if fileName == 'mo_2012_11_24_8_55_21'
      from = 1;  to = 8750;
    end
  elseif diag_sys == 2
    if fileName == 'mo_2008_10_27_15_50_49'
      from = 1;  to = 14000;
    elseif fileName == 'mo_2009_10_20_9_12_49'
      from = 1;  to = 33000;
    elseif fileName == 'mo_2009_10_21_9_17_55'
      from = 1;  to = 43000;
    elseif fileName == 'mo_2009_11_23_9_47_4'
      from = 1;  to = 23000;
    elseif fileName == 'mo_2009_11_21_9_31_4'
      from = 9000;  to = 30000;
    end
  end
  resData = initData(from : to, :);
endfunction

function resParameters = readParametersData(filePath, fileName, fileExt, initParameters)
//******************************************************************************************
// Function that reads parameters data from a archive file.                                *
// IN:  filePath - path to the archive file                                                *
//      fileName - archive file name                                                       *
//      fileExt - extention of a archive file                                              *
//      initParameters - array of initial structures of archives with information          *
//      about parameters                                                                   *
// OUT: resParameters - array of result structures of archives with readed parameters data *
//******************************************************************************************
  filePathName = filePath + '/' + fileName + fileExt;
  if ~isfile(filePathName)
    printf("[ERROR]: cannot read the data from the archive file:\n");
    printf("\t''%s''\n", filePathName);
    printf("\tThe reason: this file is not exist or it is not a file\n");
    abort;
  end
  
  archiveData = fscanfMat(filePathName);
  archiveData = getValidData(fileName, archiveData);
  count_params = size(initParameters, 1);
  resParameters = initParameters;
  for i = 1 : count_params
    resParameters(i).data = archiveData(:, initParameters(i).archIndexStart : initParameters(i).archIndexEnd);
  end
endfunction

// ======= OUT RESULTS =======
function saveSteadyPoints(path, fileName, reg, dt)
//*****************************************************
// Save the steady mode points values to a text file  *
// IN:  path - path to the file directory             *
//      fileName - name of the file                   *
//      reg - regime parameter data for saving        *
//      dt - 'dt' parameters data for saving          *
//*****************************************************
  Ncols = size(dt, 'c');  Nrows = size(dt, 'r');
  filePathName = path + '/' + fileName;
  
  if ~isdir(path)
    createdir(path);
  else
    mdelete(filePathName); // delete file with the same name as a current file
  end
  
  f = mopen(filePathName, 'w');
  
  // Write values
  for i = 1 : Nrows
    mfprintf(f, "%f ", reg(i));
    for j = 1 : Ncols
      mfprintf(f, "%f ", dt(i, j));
    end
    mfprintf(f, "\n");
  end
  mclose(f);
  
  // show info message
  printf("[INFO]: The steady mode points values was exported to the text file:\n");
  printf("\t%s\n", filePathName);
endfunction

function saveResToGraphicFiles(path_imagesOut, calcIdentif, ext_images)
//*********************************************************************************
// Save results to a graphics file.                                               *
// IN:  path_imagesOut - path to the out files                                    *
//      calcIdentif - information about current calculation (identification info) *
//      ext_images - extention for the graphics files                             *
// OUT: ---                                                                       *
//*********************************************************************************
  // go into the common directory (dir for saving all images)
  if(~isdir(path_imagesOut) & ~createdir(path_imagesOut))
    printf("[ERROR]: Cannot create the common dir for saving result images:\n\t%s\n", path_imagesOut);
    abort;
  end
  // go into the special directory (dir for current calculation)
  path_imagesOut = path_imagesOut + '/' + calcIdentif;
  if isdir(path_imagesOut)
    removedir(path_imagesOut);
  end
  if(~createdir(path_imagesOut))
    printf("[ERROR]: Cannot create the dir for saving result images: %s\n", path_imagesOut);
    abort;
  end
  
  // save image data
  printf("[INFO]: Export the results plotted images of the GTE''s oil''s characteristics to the files:\n");
  figureIDs = winsid();
  for i = 1 : length(figureIDs)
    h = scf(figureIDs(i));
    exportFileName = h.figure_name + '.' + ext_images;
    exportFileName = strsubst(exportFileName, ' ', ''); // spaces deleting
    filePathName = path_imagesOut + '/' + exportFileName;
    xs2png(h, filePathName);
    printf("\t%s\n", filePathName);
  end
endfunction

function saveResToTextFile(path, fileName, strTitle, strDateTime, powers, ..
                           Ncols, Nrows, .. 
                           par_regim_data, par_regim_name, par_oil_data, par_oil_names)
//***************************************************************************************************
// Save results to a text file.                                                                     *
// IN:  path - path to the out file                                                                 *
//      fileName - name of the out file                                                             *
//      strTitle - title string, that writes in the out file with information about current results *
//      powers - array of the polynomial powers                                                     *
//      Ncols - quantity of the columns of the oil data array (par_oil_data array)                  *
//      Nrows - quantity of the rows of the oil data array (par_oil_data array)                     *
//      par_regim_data - data of the regime parameter of the result characteristics (X-axes)        *
//      par_regim_name - name of the regime parameter of the result characteristics (X-axes)        *
//      par_oil_data - data of the oil parameters (Y-axes)                                          *
//      par_oil_names - names of the oil parameters (Y-axes)                                        *
// OUT: ---                                                                                         *
//***************************************************************************************************
  filePathName = path + '/' + fileName;
  
  if ~isdir(path)
    createdir(path);
  else
    mdelete(filePathName); // delete file with the same name as a current file
  end
  
  f = mopen(filePathName, 'w');
  mfprintf(f, "The results information:  %s\n", strTitle);
  mfprintf(f, "Date_time:  %s\n", strDateTime);
  mfprintf(f, "The polynomial powers: ");
  for i = 1 : Ncols
    mfprintf(f, "%i  ", powers(i));
  end
  
  // Write the parameters names
  mfprintf(f, "\n\t\t%s  ", par_regim_name);
  for i = 1 : Ncols
    mfprintf(f, "%s\t\t", par_oil_names(i));
  end
  mfprintf(f, "\n");
  
  // Write values
  for i = 1 : Nrows
    mfprintf(f, "%8.2f\t\t", par_regim_data(i));
    for j = 1 : Ncols
      mfprintf(f, "%5.2f\t\t\t", par_oil_data(i, j));
    end
    mfprintf(f, "\n");
  end
  mclose(f);
  
  // show info message
  printf("[INFO]: The results points values was exported to the text file:\n");
  printf("\t%s\n", filePathName);
endfunction

