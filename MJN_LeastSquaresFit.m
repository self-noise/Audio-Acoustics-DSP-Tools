%-------------------------------------------------------------------------%
% Simple least squares solver, with data fitted to arbitrary order
% Author    : Michael J. Newton
% Date      : 02/04/2018
%-------------------------------------------------------------------------%
% Updates   : 28/02/2024 (changed function name and altered output to
%             produce R-squared
%-------------------------------------------------------------------------%
% Inputs:
%   xData       : Column of x-axis points
%   yData       : Column of corresponding y-axis points
%   order       : 1 for linear, 2 for quadratic, etc
%
% Output:
%
%   LS_Soln     : Vector of size (1,N_order+1) with coeffs. of the LS fit
%   varargout{1}: Overall 'error' metric ('R-squared')
%   varargout{2}: Vector of y-values that lie along LS curve (i.e., the
%                 specific preductions of 'y' for each 'x' values in the
%                 input data), to allow easy/quick plotting of data return
%-------------------------------------------------------------------------%
function [LS_Soln,varargout] = MJN_LeastSquaresFit(xData,yData,N_order)

Npts = length(xData);

H = zeros(Npts,N_order+1);
H(:,1) = ones(Npts,1);

for nOrder = 1:N_order
    H(:,nOrder+1) = xData.^nOrder;
end

% The key bit that inverts the matrix and to produce gradient(s) and
% y-intercept for the fit
LS_Soln     = (H'*H)\H'*yData;

yData_tilde = H*LS_Soln;
R           = sum((yData-yData_tilde).^2);  % Residual sum of squares
SumSquares  = (length(yData)-1)*var(yData); % Total sum of squares

% R-squared error metric (lies between 0:1)
Rsq         = 1 - R/SumSquares;

% Optional outputs
varargout{1} = Rsq;
varargout{2} = yData_tilde;
end