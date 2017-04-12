function [fitresult, gof] = createFit(angle, x)
%CREATEFIT(ANGLE,X)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : angle
%      Y Output: x
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-Jul-2016 17:17:07


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( angle, x );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
fitresult([.1 .2 .3])


