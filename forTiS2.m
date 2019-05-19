prompt = {'Reference number of Data in Exp group:',...
    'Enter the folder of Masks for the Data:',...
    'Enter the matlab file of Data timeline:',...
    'Enter the scan rate of cyclic voltammetry (mV):',...
    'Save route:',...
    'Sampling Rate (fps)'};

dlg_title = 'Input';
num_lines = 1;

defaultans = {'A1',...
    'H:\TaS2\20190415_TiS2_Au_CH18SH\Mask\A1',...
    'H:\TaS2\20190415_TiS2_Au_CH18SH\Timer\A1_data.mat',...
    '100',...
    'H:\TaS2\20190415_TiS2_Au_CH18SH\Result',...
    '106'};

answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
expName = answer{1};
maskPath = answer{2};
loader = answer{3};
rate = str2num(answer{4});
saveRoute = answer{5};
Fs = str2num(answer{6});

% TaS2(expName, maskPath, loader, rate, saveRoute, Fs);
TiS2(expName, maskPath, loader, rate, saveRoute, Fs);