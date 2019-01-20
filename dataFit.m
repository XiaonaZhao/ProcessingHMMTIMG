function Yfit = dataFit(Y)

fo = fitoptions('method','SmoothingSpline','SmoothingParam',0.05);
ft = fittype('smoothingspline'); 
data = Y;
x_1 = (1:numel(data))';
ok = isfinite(x_1) & isfinite(data(:));
% data=(data)';
cf = fit(x_1(ok),data(ok),ft,fo);
data_fit = cf(x_1);
Yfit = data_fit;