function [tifSeq, tifNum] = ReadTifFiles(TopName, skipNum)

[tifFiles, tifPath] = uigetfile('*.tiff', 'Multiselect', 'on', TopName); 
if isequal(tifFiles, 0)
   disp('User selected Cancel')
   return;
end
tifFiles = cellstr(tifFiles);  % Care for the correct type 

if skipNum == 0
    tifNum = length(tifFiles);
else
    tifNum = floor((length(tifFiles))/(skipNum+1))+1;
end

tifSeq = cell(tifNum,1); % Line Data structure, imgSeq is column vector
for j = 1:tifNum
    tifSeq{j} = imread(fullfile(tifPath, tifFiles{j})); % read all tiff files
end
clear jj skipNum
end