// Additional functions

function [x, y] = sortByX(x, y)
//*********************************************************************
// Sorting by "x" array. Others arrays, that stores in "y" variable,  *
// will be sorted in accordance with with "x" array.                  *
//*********************************************************************

  n = length(x);
  col = size(y, 'c');
  for j = 1 : n
    for i = 2 : n
      if (x(i) < x(i - 1))
        
        tempX = x(i);
        x(i) = x(i - 1);
        x(i - 1) = tempX;
        
        for k = 1 : col
          tempY(k) = y(i, k);
          y(i, k) = y(i - 1, k);
          y(i - 1, k) = tempY(k);
        end
      
      end
    end
  end  
endfunction

// Mathematical functions
function yu = interExtraPolation(x, y, xu)
//********************************************************
// The function of the linear inter- and extrapolations  *
//********************************************************

  n = length(x);
  i = 2;
  while (xu > x(i))
    i = i + 1;
    if (i > n)
      x1 = x(n - 1);
      x2 = x(n);
      y1 = y(n - 1);
      y2 = y(n);
      yu = (xu - x1) * (y2 - y1) / (x2 - x1) + y1;
      return;
    end
  end
  if (xu == x(i))
    yu = y(i);
    return;
  elseif (xu < x(i))
    x1 = x(i - 1);
    x2 = x(i);
    y1 = y(i - 1);
    y2 = y(i);
    yu = (xu - x1) * (y2 - y1) / (x2 - x1) + y1;
  end
endfunction


// Approximation
function coef = coeffs_trend_n(x, y, n)
//*****************************************************************
// The function, that calculate coefficients of n-degree trend    *
// In:  x, y - the arrays of points, that approximate             *
//      n - the polynomial power value                            *
// Out: the trend coefficients: a, b, c, d ... - coef(1 2 3 4...) *
//*****************************************************************

  M = 2 * n + 1;
  
  for j = 1 : n + 1
    for i = 1 : n + 1
      K(j, i) = sum(x .^ (M - i - j + 1));
    end
  end
  K(n + 1, n + 1) = length(x);

  for i = 1 : n + 1
    S(i) = sum(x .^ (n + 1 - i) .* y);
  end
  
  coef = inv(K) * S;
endfunction

function y_apr = approximation(x_init, y_init, x_apr, powers, Nrows_apr, Ncols_apr)
//****************************************************************************************************
// Function, that perform approximation of distributions points by polynomial line with any power.   *
// The main features are:                                                                            *
// - perform the calculation of approximate values for "x_apr" array values;                         *
// - if arrays "x_init", "y_init" and "x_apr" is more than 1-dimension arrays (need to               *
//   approximate more than 1 line), every line can be approximated by the different powers,          *
//   stores in the "powers" array variable.                                                          *
// In:  x_init - array of the "x" initial points values                                              *
//      y_init - array of the "y" initial points values                                              *
//      x_apr - array of the "x" approximation line points values                                    *
//      powers - array of the polynomial powers (for every approximation line)                       *
// Out: y_apr - array of the "y" polynomial line points values                                       *
//****************************************************************************************************
  y_apr(Nrows_apr, Ncols_apr) = 0;
  for i = 1 : Ncols_apr
    coefs = coeffs_trend_n(x_init, y_init(:, i), powers(i)); // define the polynomial coefficients for the current line
    for j = 1 : powers(i) + 1
      y_apr(:, i) = y_apr(:, i) + coefs(j) * x_apr .^ (powers(i) + 1 - j);
    end
  end
endfunction

// Data & Time
function date_time_str = getDateTimeString()
//********************************************************************
// Preparing the string in the format: "year.month.day_hour:min:sec" *
//********************************************************************
  date_time_int = int(datevec(now()));
  date_time_str = sprintf("%i.%i.%i_%i:%i:%i", date_time_int(1), date_time_int(2),..
                                               date_time_int(3), date_time_int(4),..
                                               date_time_int(5), date_time_int(6));
endfunction

