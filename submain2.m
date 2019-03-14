%% fullfill the expTab (ROImean, inside & outside mean intensity)
tic

matPath = 'G:\Summary\MoS2';

load 'HMMT_tif_mask_FileName.mat'
[foldernumber, col] = size(expTab);

if foldernumber ~= size(maskTab, 1)
    return
end

addCell = cell(foldernumber, 8);
expTab = [expTab addCell];

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for ii = 1:1:foldernumber
    maskFile = convertStringsToChars(strcat(maskTab{ii, 2}, maskTab{ii, 3}));
    mask = ~imread(maskFile);
    tifFile = convertStringsToChars(strcat(expTab{ii, 4}, expTab{ii, 5}));
    tif = ReadTifROI(tifFile, mask);
    
    [expTab{ii, col+1}, expTab{ii, col+2}] = calculateDonutEdge(tif{expTab{ii, 2}, 3}, mask);
    expTab{ii, col+3} = expTab{ii, col+2} - expTab{ii, col+1};
    [expTab{ii, col+4}, expTab{ii, col+5}] = calculateDonutEdge(tif{expTab{ii, 3}, 3}, mask);
    expTab{ii, col+6} = expTab{ii, col+5} - expTab{ii, col+4};
    expTab{ii, col+7} = expTab{ii, col+6} - expTab{ii, col+3};
    
    expTab{ii, col+8} = tif(:, 4);
    
    cellpath = [matPath '\' '2xHMMT_' mat2str(ii) '.mat'];
    save(cellpath, 'tif', '-v7.3');
    clear tif
    
    save('HMMT_expTab_FileName.mat', 'expTab');
    
    if foldernumber - ii <= 1
        waitbar(ii/foldernumber, hwait, 'Test to be accompleted');
        pause(0.05);
    else
        PerStr = round(ii/foldernumber*100);
        str = ['Running', num2str(PerStr), '%'];
        waitbar(ii/foldernumber, hwait, str);
        pause(0.05);
    end
    
end
close(hwait);

MailToMe('nona1588@outlook.com');

toc
disp(['Running time: ',num2str(toc)]);

warndlg('Test passed', 'Warning')

%% detrend ROImean
addCell = cell(foldernumber, 3);
expTab = [expTab addCell];
for ii = 1:1:foldernumber
    frames = size(expTab{ii, col+8}, 1);
    x = (1:1:frames)';
    
    if ii == 3 || 4
        Fs = 200;
    else
        Fs = 100;
    end
    
    rawCurve = cell2mat(expTab{ii, col+8});
    filterCurve = lowp(rawCurve, 4, 30, 0.1, 20, Fs);
    
    expTab{ii, col+9} = driftBaseline(x, filterCurve);
    pause(1);
    
    expTab{ii, col+10} = expTab{ii, col+9}(expTab{ii, 2});
    expTab{ii, col+11} = expTab{ii, col+9}(expTab{ii, 3});
    
    close all
end

save('HMMT_exp_FileInfo.mat', 'expTab');

%% manage 'expTab'
fields = {'number',... % 1
    'Darked',... % 2
    'Reduced',... % 3
    'disc',... % 4
    'index',... % 5
    'RingDarkedInside',... % 6
    'RingDarkedOutside',... % 7
    'RingDarkedDelta',... % 8
    'RingReducedInside',... % 9
    'RingReducedOutside',... % 10
    'RingReducedDelta',... %11
    'RingReDvsDaD',... %12 = 11-8
    'ROImean',... % 13
    'ROImeanDetrended',... % 14
    'ROImeanDarked',... % 15
    'ROImeanReduced' %16
    };
expStruct = cell2struct(expTab, fields, 2);

save('HMMT_exp_FileInfo.mat', 'expStruct', '-append');

%%  intensity to concentrate

Fs = 100;
filterCurve = lowp(rawCurve, 4, 30, 0.1, 20, Fs);
frames = size(filterCurve, 1);
x = (1:1:frames)';
treatedCurve = driftBaseline(x, filterCurve);

standard = -treatedCurve(234+1);
BalphaOxi = standard/10; % mM^(-1)

foldernumber = size(expTab, 1);
redConcentrate = zeros(foldernumber, 1);

for ii = 1:1:foldernumber
    if expStruct(ii).ROImeanReduced > standard
        redConcentrate(ii) = 10;
    elseif expStruct(ii).ROImeanReduced < standard && expStruct(ii).ROImeanReduced > 0
        redConcentrate(ii) = expStruct(ii).ROImeanReduced/BalphaOxi;
    else
        disp(['***' 'Please check data in expStruct' ii 'at ROImeanReduced***']);
    end
end

