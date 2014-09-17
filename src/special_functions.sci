// SPECIAL FUNCTIONS

// The power (N) defining
function N = p2ToPower(p2, count, p2_characs, N_characs)
//*****************************************************************************************
// Define the GTE's power values (Ngte) with yhe help of the characteristics Ngte = f(p2) *
// IN:  p2 - the array with the p2 parameter values                                       *
//      count - the "p2" array size                                                       *
//      p2_characs - the x-values (p2) of the characteristics for defining powers         *
//      N_characs - the x-values (Ngte) of the characteristics for defining powers        *
// OUT: N - the array with the result Ngte parameter values                               *
//*****************************************************************************************
  for i = 1 : count
    N(i) = interExtraPolation(p2_characs, N_characs, p2(i));
  end
endfunction

function N = ngvToPower(ngv, N_nom, ngv_nom)
//*********************************************************************************************
// Define the power on the outlet reduction shaft (Nred) with the help of the reducer outlet  *
// shaft rotation speed (ngv) and the theoretical cubic dependence Nred and ngv.              *
// IN:  ngv - array of the reducer outlet shaft rotation speeds values                        *
//      N_nom - the value of Ngte parameter on the nominal works regime                       *
//      ngv_nom - the value of ngv on the nominal works regime                                *
// OUT: N - array of the powers on the reducer outlet shaft                                   *
//*********************************************************************************************
  power_N_ngv = 3;
  a = (N_nom / ngv_nom ^ power_N_ngv);
  N = a * ngv ^ power_N_ngv;
endfunction

function [invalidCols, invalidRows] = calcInvalidValuePos(indexes, count_rows)
//*****************************************************************************************
// Calculate the positions (rows, cols) of the parameters arrays rows with invalid data   *
// IN:  indexes - array of indexes with invalid values of parameters                      *
//      count_rows - the rows quantity in parameters arrays                               *
// OUT: invalidRows - the array of the parameters arrays rows numbers with invalid values *
//      invalidCols - the array of the parameters arrays cols numbers with invalid values *
//*****************************************************************************************
  invalidCols = ceil(indexes / count_rows);
  invalidRows = indexes - int(indexes / count_rows) * count_rows;
  lastValuesInCols = find(invalidRows == 0);
  if length(lastValuesInCols)
    invalidRows(lastValuesInCols) = count_rows;
  end
endfunction

function plotInvalidArchive( reg, tm, reg_steady, tm_steady, arrayNumber_steady, ..
                             cols_invalid_dt, rows_invalid, rows_invalid_u, ..
                             index_in, reg_name, t_name, ..
                             str_archiveNumberName, colors )
//*************************************************************************************************************
// Function for plotting archive parameters with the invalid steady mode points.                              *
// IN:  reg - the regime parameter time series                                                                *
//      tm - the temperature parameters time series                                                           *
//      reg_steady - the regime parameter steady mode points                                                  *
//      tm_steady - the temperature parameter steady mode points                                              *
//      arrayNumber_steady - the steady mode points of array numbers                                          *
//      cols_invalid_dt - the number(-s) of columns with invalid point(-s) value(-s) (if use 'dt' parameter)  *
//      rows_invalid - the number(-s) of rows with invalid point(-s) value(-s)                                *
//      rows_invalid_u - the unique values from the 'rows_invalid' arrays                                     *
//      index_in - the number of the index temperature at the entry                                           *
//      reg_name - the name of the regime parameter                                                           *
//      t_name - the array of the temperatures names                                                          *
//      str_archiveNumberName - the string of union information about number and name archive                 *
//      colors - the colors indexes arrays                                                                    *
// OUT: ---                                                                                                   *
//*************************************************************************************************************
  // define the valid steady mode points
  reg_steady_valid = reg_steady;  reg_steady_valid(rows_invalid_u, :) = [];
  tm_steady_valid = tm_steady;  tm_steady_valid(rows_invalid_u, :) = [];
  arrayNumber_steady_valid = arrayNumber_steady;  arrayNumber_steady_valid(rows_invalid_u, :) = [];
  
  // Initial data for plotting
  x_time = [0 : length(reg) - 1]; // X-value for plotting
  kRegScale = 0.1; // scaled coefficient for plotting 'reg' parameter in same window as 'tm' parameters
  type_validPoints = -3; // type of the marker for plotting the valid points
  type_invalidPoints = -9; // type of the marker for plotting the invalid points
  legend_str = [];
  cols_invalid_t = cols_invalid_dt + 1; // conversion columns number from 'dtm' parameter to 'tm'
  cols_invalid_t_u = unique(cols_invalid_t); // numbers of the parameters with invalid steady mode points
  count_invalidParams = length(cols_invalid_t_u); // quantity of the parameters with invalid steady mode points

  // plot window settings
  windowNumber = max(winsid()) + 1;
  strTitle = 'Invalid points, ' + str_archiveNumberName;
  hPlot = scf(windowNumber); xgrid;
  hPlot.figure_name = strTitle;
  title(strTitle, 'fontsize', 4);
  
  // plot parameters time series graphs
  plot2d(x_time, reg * kRegScale, colors(1));  e = gce(); e.children.thickness = 2;
  legend_str = [legend_str; reg_name];
  plot2d(x_time, tm(:, index_in), colors(2));  e = gce(); e.children.thickness = 2;
  legend_str = [legend_str; t_name(index_in)];
  for i = 1 : count_invalidParams
    plot2d(x_time, tm(:, cols_invalid_t_u(i)), colors(i + 2));  e = gce(); e.children.thickness = 2;
    legend_str = [legend_str; t_name(cols_invalid_t_u(i))];
  end
  legend(hPlot, legend_str, 1);

  // plot invalid steady mode points
  plot2d(arrayNumber_steady(rows_invalid_u), reg_steady(rows_invalid_u) * kRegScale, type_invalidPoints);
  plot2d(arrayNumber_steady(rows_invalid_u), tm_steady(rows_invalid_u, index_in), type_invalidPoints);
  for i = 1 : count_invalidPoints
    plot2d(arrayNumber_steady(rows_invalid(i)), tm_steady(rows_invalid(i), cols_invalid_t(i)), type_invalidPoints);
  end

  // plot valid steady mode points
  if length(reg_steady_valid) ~= 0
    plot2d(arrayNumber_steady_valid, reg_steady_valid * kRegScale, type_validPoints);
    plot2d(arrayNumber_steady_valid, tm_steady_valid(:, index_in), type_validPoints);
    for i = 1 : count_invalidParams
      plot2d(arrayNumber_steady_valid, tm_steady_valid(:, cols_invalid_t_u(i)), type_validPoints);
    end
  else
    printf("[INFO]: None steady mode, founded in current archive, is valid\n");
  end
endfunction

