// Special functions
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

