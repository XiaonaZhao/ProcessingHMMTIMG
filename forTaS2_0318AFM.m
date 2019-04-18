prompt = {'Reference number of Data in Exp group:',...
    'Enter the folder of Masks for the Data:',...
    'Enter the matlab file of Data timeline:',...
    'Enter the scan rate of cyclic voltammetry (mV):',...
    'Save route:'};

dlg_title = 'Input';
num_lines = 1;

defaultans = {'A1',...
    'F:\TaS2\20190324_TaS2_0318_ITO\MaskA1',...
    'G:\TaS2\20190318_TaS2_ITO\matlab_data\B2.mat',...
    '100',...
    'F:\TaS2\20190324_TaS2_0318_ITO\result'};

answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
expName = answer{1};
maskPath = answer{2};
loader = answer{3};
rate = str2num(answer{4});
saveRoute = answer{5};

% if DC
TaS2_DC(expName, maskPath, saveRoute);
% if CV
% TaS2(expName, maskPath, loader, rate, saveRoute);

%%
Fs = 106;
y = lowp(a(:, 2), 2, 43, 0.1, 20, Fs);
plot(x, y)
%%
load('F:\TaS2\20190324_TaS2_0318_ITO\result\expTab_TaS2_0324.mat')

for m = 1:size(expTab_TaS2, 1)
    expName = expTab_TaS2{m, 2};
    tifPath = expTab_TaS2{m, 3};
    maskPath = expTab_TaS2{m, 4};
    saveRoute = expTab_TaS2{m, 5};
    TaS2_DC(expName, tifPath, maskPath, saveRoute);
end
