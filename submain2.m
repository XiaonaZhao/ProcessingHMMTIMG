matPath = 'J:\Summary\MoS2';

load 'J:\Summary\MoS2\HMMT_tif_mask_FileName.mat'
[foldernumber, col] = size(expTab);

if foldernumber ~= size(maskTab, 1)
    return
end

addCell = cell(foldernumber, 9);
expTab = [expTab addCell];
tifFile = expTab(:, 4);
hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for ii = 1:1:foldernumber
    mask = ~imread(maskTab{ii, 2}); % Note: maskTab{:, 2} is string not cell.
    tif = ReadTifROI(tifFile{ii}, mask);
    
    expTab{ii, col+8} = tif(:, 4);
    
    cellpath = [matPath '\' '2xHMMT_' mat2str(ii) '.mat'];
    save(cellpath, 'tif', '-v7.3');
    clear tif
    
    [expTab{ii, col+1}, expTab{ii, col+2}] = calculateDonutEdge(expTab{ii, 2}, mask);
    expTab{ii, col+3} = expTab{ii, col+2} - expTab{ii, col+1};
    [expTab{ii, col+4}, expTab{ii, col+5}] = calculateDonutEdge(expTab{ii, 3}, mask);
    expTab{ii, col+6} = expTab{ii, col+5} - expTab{ii, col+4};
    expTab{ii, col+7} = expTab{ii, col+6} - expTab{ii, col+3};
    
    if foldernumber - ii <= 5
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
warndlg('Test passed', 'Warning')
