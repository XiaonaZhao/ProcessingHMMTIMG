Fs = 100;
rate = 50;

load('E:\Pollutions\MoS2_20200718_ITO_pNP\_Timer\B2_data.mat')
begin = triggerTime_AC(data, t, Fs);

tifFile = 'E:\Pollutions\MoS2_20200718_ITO_pNP\B2_pNP_0_-0-8V_50mVpS_2c_Pike100fps';
tifDir = dir(fullfile(tifFile, '*.tiff'));

potential = potentialLine(rate, Fs, 0, -0.8); 


validDir = tifDir(begin.frame:(begin.frame+length(potential)));

[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));

tifPage = zeros(length(potential), length(sROI));
tif0 = double(imread(fullfile(tifFile, validDir(1).name)));
for ii = 1:length(potential)
    tif1 = double(imread(fullfile(tifFile, validDir(ii).name)));
    tif2 = double(imread(fullfile(tifFile, validDir(ii+1).name)));
    tif = (tif2 - tif1)./tif0;
    for jj = 1:length(sROI)
        [row, col] = ImageJroiLocation(sROI{jj});
        temp = tif(row(1):row(2), col(1):col(2));
        tifPage(ii, jj) = mean(temp(:));
    end
end
%%
curve = lowp(tifPage(:, 3), 1, 5, 0.1, 20, 100);

% figure('color', 'w');
plot(potential, curve)