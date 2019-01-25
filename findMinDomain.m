function centralBW = findMinDomain(BW)

imLabel = bwlabel(BW);% Mark each connected domain
% imLabel = logical(BW);% it doesn't work
stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
area = cat(1, stats.Area);
index = find(area == min(area));% Find the index of the smallest connected domain
centralBW = ismember(imLabel, index);