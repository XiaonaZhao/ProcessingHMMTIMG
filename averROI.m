function Intensity = averROI(roiSeq, roiNum, nonzeroBW)

Intensity = zeros(roiNum, 1);
for j = 1:roiNum
    temp = roiSeq{j}; % temp is 640x480
    Intensity(j) = sum(temp(:))/nonzeroBW; % Intensity is column vector
end
clear j
end