function [tifFolder, tifNames] = ReadTifFileNames(tifFile)

tifFolder = fullfile(tifFile);
dirOutput = dir(fullfile(tifFolder, '*.tiff'));
dirFileTif = sortObj(dirOutput);
tifNames = {dirFileTif.name}';
end