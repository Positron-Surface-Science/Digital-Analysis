function [fitresult, gof] = createFit(x, undo2)

[xData, yData] = prepareCurveData( x, undo2 );

% Set up fittype and options.
ft = fittype( '-2.62E-6*x+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.544727506301423;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end
