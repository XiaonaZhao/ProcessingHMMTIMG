% the intensity vs. Potential at a high frequency
load('E:\MoS2\MoS2_0919_0802\_Timer\A3_data.mat')
begin = triggerTime_PEIM(data, t, 100);
%%
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
sROI = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));

%%
load('E:\MoS2\MoS2_0919_0802\_Result\A_sequance\A3.mat')
%%
[Value.tifPath, Value.tifName] = ReadTifFileNames(Value.tifFile);
tif1 = im2double(imread(fullfile(Value.tifPath, Value.tifName{1})));
L = 1024;

average = zeros(L, 1);
for jj = begin.frame:(begin.frame+L-1)
    tif = im2double(imread(fullfile(Value.tifPath, Value.tifName{jj}))) - tif1;
    
    [row, col] = ImageJroiLocation(sROI{1});
    tif0 = tif(row(1):row(2), col(1):col(2));
    average(jj-begin.frame+1, 1) = mean(tif0(:));
    
end
%%
voltage = data(begin.CS1:begin.CS1+1023*100, 2);
voltage = voltage(1:100:end);
%%
h = voltage-0.8*ones(1024, 1);
%%
t = (0:1023)'/100;
plot(t, 100*average+0.3*ones(1024, 1), t, h/2)
ylim([-0.3 0.3])