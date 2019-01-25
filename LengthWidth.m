% Projection aspect ratio of irregular graphics
L = zeros(row,1);
B = zeros(row,1);
LvsB = zeros(row,1);
for n = 1:row
    img0 = imread(fullfile(fileFolder, fileNames{n}));
    imLabel = bwlabel(img0);% Mark each connected domain
    stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
    area = cat(1, stats.Area);
    index = find(area == min(area));% Find the index of the smallest connected domain
    img = ismember(imLabel, index);
    k = find(sum(img));
    L(n) = k(end) - k(1) + 1;
    j = find(sum(img, 2));
    B(n) = j(end) - j(1) + 1;
    LvsB(n) = L(n)/B(n);
end