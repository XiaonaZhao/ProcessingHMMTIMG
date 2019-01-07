function imgSubtractBg = BgdRemoval(imgSeq, imgNum,BgSeq, BgNum)

sum = zeros(480, 640);
for j = 1:BgNum
    sum = double(BgSeq{j}) + sum;
end
averBg = sum./double(BgNum);
clear j BgSeq BgNum

imgSubtractBg = cell(imgNum,1);
for j = 1:imgNum
    imgSubtractBg{j} = double(imgSeq{j}) - averBg; 
end
clear j imgSeq imgNum
end