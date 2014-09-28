// SPECIAL FUNCTIONS

function [reg_all, dtm_all] = calcSteadyPoints()
//*****************************************************************
// Function for calculating the steady mode points values         *
// IN:  all initial data is the global variables                  *
// OUT: reg_all - the regime parameter steady mode points values  *
//      dtm_all - the 'dt' parameters steady mode points values   *
//*****************************************************************
  // Arrays for storing steady mode points values
  reg_all = [];
  dtm_all = [];
  
  // MAIN CYCLE
  for fileIndex = 1 : size(filesArchive, 'r')
    // Deleting results that was obtained in the previous main cycle iteration
    clear reg_steady; clear tm_steady; clear dtm_steady; clear arrayNumber_steady;
  
    //  READING an archive file
    params = readParametersData(path_archives, sep, filesArchive(fileIndex), '.' + ext_archive, params);

    // Getting the parameters arrays with reduction to a normal atmospheric condition
    // This parameters is take out of the if-else, because name of its index variable for the both diagnostic systems is equal
    reg = params(index_reg).data; // the regime parameter values
    tm = params(index_tm).data;
    if diag_sys == GTE_OIL
      t0 = mean(params(index_t0).data, 'c');
      alpha = sqrt(288 ./ (t0 + 273)); // alpha coefficient for reductions parameters to normal atmospheric conditions
      Gt = params(index_Gt).data;
      reg = reg * 10.2;  // there are p2 parameter values with conversion its from MPa to kg/cm2
      n2 = params(index_n2).data .* alpha; // reduction to normal conditions
    end
    
    //  Splitting parameters arrays to sectors with defining the average values and strange
    steadyIndex = 0;
    // kk - is the coefficient for correct definition the iterations quantity and correct shifting buffer with length "sectorShift" 
    //      along the full length of archive
    kk = ceil(sectorLength / sectorShift);
    arrSize = length(reg);
    for j = kk : int(arrSize / sectorShift)
      from = sectorShift * (j - 1) + 1;
      to = sectorShift * j;
      arrayNumber = j * sectorShift;
      
      // split arrays, calc strange or average values of parameters and define steady mode
      isSteadyMode = %F;
      if diag_sys == GTE_OIL
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
        if diag_sys == GTE_OIL
          reg_steady(steadyIndex) = median(reg(from : to));
        else
          reg_steady(steadyIndex) = reg_avrg; // already calculated
        end
        //--------------------------------------------------------------------------------------------------------
        // calculate the 'tm' values in steady mode of work
        forecastTo = to + forecastInterval; // the argument-value for forecasting
        xModel = [arrayNumber - modelLength + 1 : arrayNumber]';
        yModel = tm(from : to, :);
        tm_steady(steadyIndex, :) = linearForecastValues(xModel, yModel, forecastTo)';
        
        // old version
//        for t = 1 : count_tmParams
//          tm_steady(steadyIndex, t) = median(tm(from : to, t));
//        end
        //--------------------------------------------------------------------------------------------------------
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
        abort;
      end
    end

    //  Define temperature drop of oil
    for t = 1 : count_dtmParams
      dtm_steady(:, t) = tm_steady(:, t + 1) - tm_steady(:, index_in);
    end

    // Check existence the "bad", invalid points, that is far from others points
    ind_invalidValues = find(dtm_steady > Udtm_valid_max | dtm_steady < Udtm_valid_min);
    count_invalidPoints = length(ind_invalidValues);
    if count_invalidPoints
      [cols_invalid_dt, rows_invalid] = calcInvalidValuePos(ind_invalidValues, size(reg_steady, 'r'));
      rows_invalid_u = unique(rows_invalid);
      str_archiveNumberName = 'archive #' + string(fileIndex) + ': ' + filesArchive(fileIndex)';

      printf("[ERROR]: There was found %i invalid steady mode point(s) in the %s\n", count_invalidPoints, str_archiveNumberName);
      for i = 1 : count_invalidPoints
        printf("\tpoint #%i: parameter = ''%s'', number = %i\n", i, dt_names(cols_invalid_dt(i)), rows_invalid(i));
      end
      
      printf("Continue with deleting invalid points? (1 - yes, 2 - no)\n");
      key = scanf("%i");
      if key == 1
        // delete rows with invalid point(-s) for getting arrays with only valid points
        reg_steady(rows_invalid_u, :) = [];
        dtm_steady(rows_invalid_u, :) = [];
      else
        plotInvalidArchive( reg, tm, reg_steady, tm_steady, arrayNumber_steady, ..
                            cols_invalid_dt, rows_invalid, rows_invalid_u, ..
                            index_in, params(index_reg).name, t_name, ..
                            str_archiveNumberName, colors );
        printf("--------------------------------------------------------------------------------------------------------\n\n");
        abort;
      end
    end

    // Save the values of the steady mode points over all archives for further processing out of the main cycle
    reg_all = [reg_all; reg_steady];
    dtm_all = [dtm_all; dtm_steady];
    printf("[INFO]: Archive #%i: ""%s"", points quantity = %i\n", fileIndex, filesArchive(fileIndex), steadyIndex);
  end
endfunction

function forc_y = linearForecastValues(xArr, yArr, forc_x)
//**********************************************************************************************
// Calculation the forecast value Y for given model and argument value X (linear forecasting)  *
// IN:  xArr - array values of model                                                           *
//      yArr - array values of model                                                           *
//      forc_x - the x-argument value                                                          *
// OUT: forc_y - the forecasted value                                                          *
//**********************************************************************************************
  Ncols = size(yArr, 'c');  power = 1;
  forc_y(Ncols) = 0;
  for i = 1 : Ncols
    coefs = coeffs_trend_n(xArr, yArr(:, i), power);
    for j = 1 : power + 1
      forc_y(i) = forc_y(i) + coefs(j) * forc_x .^ (power + 1 - j);
    end
  end
endfunction

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

// GRAPHS PLOT
function plotResults(x_points, y_points, x_polyn, y_polyn, ..
                      count_pars, plotInSameWin, x_names, y_name, strDateTime)
//************************************************************
// 
//************************************************************
  // define the size of graphs square matrix
  if plotInSameWin
    size_graphsSquareMatrix = 2;
  else
    size_graphsSquareMatrix = 1;
  end
  graphsOnWin = size_graphsSquareMatrix * size_graphsSquareMatrix;
  
  // plotting
  str_y_name = ' = f(' + y_name + ')';
  oneWin_parNames = '';    sep_win_par_names = ', ';
  number_winFirst = max(winsid()) + 1;
  it_par = 1;   // iterator of plotted parameters
  it_win = 1;   // iterator of plot windows
  it_graph = 1; // iterator of the graphs on a one plot window
  while (it_par <= count_pars)
    strTitle = x_names(it_par) + str_y_name;
    oneWin_parNames = oneWin_parNames + x_names(it_par) + sep_win_par_names; // for the windows names
    hWin = scf(number_winFirst + it_win); 
    subplot(size_graphsSquareMatrix, size_graphsSquareMatrix, it_graph);
    plot2d(x_points, y_points(:, it_par), -9); xgrid;
    title(strTitle + ',  ' + strDateTime, 'fontsize', 4);

    //  Plot the approximate line
    plot2d(x_polyn, y_polyn(:, it_par), 15); e = gce(); e.children.thickness = 2;
    // Plot the approximate lines points
    plot2d(x_polyn, y_polyn(:, it_par), -14);

    if (it_graph == graphsOnWin) | (it_par == count_pars)
      // forming the windows names
      oneWin_parNames = part(oneWin_parNames, [1 : length(oneWin_parNames) - length(sep_win_par_names)]);
      hWin.figure_name = '[' + oneWin_parNames + ']' + str_y_name;
      oneWin_parNames = '';
      hWin.figure_size = [1000 700];
      hWin.figure_position = [50 50];
      
      // increase the iterators values
      it_win = it_win + 1;
      it_graph = 1;
    else
      it_graph = it_graph + 1;
    end
    it_par = it_par + 1;
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
  cols_invalid_t = cols_invalid_dt + 1; // conversion columns number from 'dt' parameter to 'tm'
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

