function [x, y] = oCenter(img)
% https://blog.csdn.net/m0_37816922/article/details/85047333
[m,n] = size(img);
sumImg = sum(img(:));
x = sum(img)*(1:n)'/sumImg;
y = (1:m)*sum(img,2)/sumImg;

end