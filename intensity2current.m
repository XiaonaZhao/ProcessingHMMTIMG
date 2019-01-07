function Current = intensity2current(Intensity, imgNum)
% The reaction on the electrode surface: H2O + e- = H + OH-

n = -1;    % Number of electrons transferred in the cell;
F = 9.65e4;   % Faraday constant;
D = 1e-9; % Actually, Do & Dr are ~ (5 - 8)e-9

AlphaOxi = 7.34;
AlphaRe = 7.34./2;
Constant = (-1)*(AlphaOxi-AlphaRe)/(D^(0.5))/(n*F*(pi^(0.5)));
% frameRate = 106; % PIKE
frameRate = 100; % HAMAMATSU
Time = (1:imgNum)./frameRate;
DeltaT = ones(imgNum,1)./frameRate;

V = ones(imgNum,1)*Time;
H = (ones(imgNum,1)*[Time(3:end) 0 0])';
Sub = sqrt(H-V);

fo = fitoptions('method','SmoothingSpline','SmoothingParam',0.05);
ft = fittype('smoothingspline'); 
data = Intensity;
x_1 = (1:numel(data))';
ok = isfinite(x_1) & isfinite(data(:));
% data=(data)';
cf = fit(x_1(ok),data(ok),ft,fo);
data_fit = cf(x_1);
Intensity = data_fit;

Current = zeros(imgNum-1,1);
Sum_new = zeros(imgNum-2,1);
for i = 3 : imgNum
    Sum_new(i) = sum(real(Current(1:i-2)'./Sub(i-2,1:i-2))).*(DeltaT(i));
    Current(i-1) = (Intensity(i-1)./Constant - Sum_new(i))./abs(DeltaT(i-1)).*((Time(i)-Time(i-1)).^(0.5));
end
Current = Current*(-1);
end

