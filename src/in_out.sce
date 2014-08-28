// Functions for performing the in-, out- operations

// ======= DATA READING =======
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
  archiveData = fscanfMat(filePath + fileName + fileExt);
  archiveData = getValidData(fileName, archiveData);
  count_params = size(initParameters, 1);
  resParameters = initParameters;
  for i = 1 : count_params
    resParameters(i).data = archiveData(:, initParameters(i).archIndexStart : initParameters(i).archIndexEnd);
  end
endfunction

// ======= OUT RESULTS =======
// Data & Time
function saveResultsToTextFile(path, fileName, strTitle, powers, ..
                               Ncols, Nrows, .. 
                               par_regim_data, par_regim_name, par_oil_data, par_oil_names)
//*******************************
//  Save results to a text file *
//*******************************
  filePathName = path + '/' + fileName;
  
  if ~isdir(path)
    createdir(path);
  else
    mdelete(filePathName); // delete file with the same name as a current file
  end
  
  f = mopen(filePathName, 'w');
  mfprintf(f, strTitle + "powers = ");
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

