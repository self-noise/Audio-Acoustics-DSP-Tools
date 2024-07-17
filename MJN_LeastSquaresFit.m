%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%   [LS_Soln,varargout] = MJN_LeastSquaresFit(xData,yData,order)
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             April 2018
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%   Performs a simple least-squared regression on a 1-dimensional dataset (x- and y-values 
%   given by the user), with a customisable order of fit (linear, quadratic, etc)
%---------------------------------------------------------------------------------------------%
% INPUTS:
%   xData       : Column of x-axis points
%   yData       : Column of corresponding y-axis points
%   order       : 1 for linear, 2 for quadratic, etc
%
% OUTPUTS:
%   LS_Soln     : Vector of size (1,N_order+1) with coeffs. of the LS fit
%   varargout{1}: Overall 'error' metric ('R-squared')
%   varargout{2}: Vector of y-values that lie along LS curve (i.e., the
%                 specific preductions of 'y' for each 'x' values in the
%                 input data), to allow easy/quick plotting of data return
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   NOTE 1: T
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   TBC
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%   2024-02-28:     Changed function name and altered output to produce R-squared value as first
%                   optional output argument
%   2018-04-02:     Created this function as a quick way to do linear data fit 'by hand'
%---------------------------------------------------------------------------------------------%
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