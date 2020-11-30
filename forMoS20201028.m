%% in this data, the sROI in expTab was for the Au
L = 1024;
Fs = 100; % sample rate

for jj = 1:11
    tifPath = expTab(jj).tifPath;
    tifDir = dir(fullfile(tifPath, '*.tiff'));
    tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));
    
    Au = zeros(length(expTab(jj).roi1), 1);
    
    for ii = 1:length(expTab(jj).roi1)
        temp = double(imread(fullfile(tifPath, tifDir(ii).name))) - tif0;
        [row, col] = ImageJroiLocation(expTab(jj).sROI);
        temp2 = temp(row(1):row(2), col(1):col(2));
        Au(ii, 1) = mean(temp2(:));
    end
    
    expTab(jj).Au = Au;
end


%%
% Amp = zeros(size(expTab, 1), 1);
% Ph = zeros(size(expTab, 1), 1);
for ii = 1:size(expTab, 1)
    begin = expTab(ii).begin;
    roi1 = expTab(ii).roi1;
%     Au = expTab(ii).Au;
    data = expTab(ii).data;
    voltage = data(:, 2);
    [expTab(ii).Amp, expTab(ii).Ph] = line(begin, roi1, voltage);
end

function [Amp, Ph] = line(begin, roi1, voltage)
L = 1024;
Fs = 100; % sample rate
voltage = voltage(begin.CS1:(begin.CS1 + L*Fs - 1)); % changed

X = roi1(begin.frame:(begin.frame + L -1));

Y = fft(X);
P2 = abs(Y/L);
Amp = max(2*P2(3:end-1));
Ph = freqrespmeasure(voltage, X);
end