function tif = ReadTifROI(fileFile, mask)
% Something wrong and I don't know how to modify this file.

[tifFolder, tifNames] = ReadTifFileNames(fileFile);
row = size(tifNames, 1);

tif = cell(row, 5);
tif0 = double(imread(fullfile(tifFolder, tifNames{1})));
for n = 1: row
    tif{n, 2} = n;
    tif{n, 3} = (double(imread(fullfile(tifFolder, tifNames{n}))) - tif0).*mask;
    tif{n, 4} = ROImean(tif{n, 2}, mask);
end