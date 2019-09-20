function [b1,reT] =binary_otus( f1 )
% otus法二值化
[r c]=size(f1);
gray_level=256;
resultT=0;
hist=zeros(gray_level,1);
for i=1:r
    for j=1:c
        hist(f1(i,j)+1)=hist(f1(i,j)+1)+1;
    end
end
hist=hist/(r*c);
vmax=0;
for tt=1:gray_level
    T=tt-1;
    w0=0;w1=0;u0=0;u1=0;var=0;%重要的地方，在每次循环遍历的时候。要把涉及到的变量清零。
    for i=1:T+1
        w0=w0+hist(i);
        u0=u0+(i-1)*hist(i);
    end
    u0=u0/w0;
    w1=1-w0;
    for j=T+2:gray_level
        u1=u1+(j-1)*hist(j);
    end
    u1=u1/w1;       
    var=w0*w1*(u0-u1)^2;
%     v(tt)=var;
    if var>vmax
        vmax=var;
        resultT=T;
    end    
end
b1=im2bw(mat2gray(f1),resultT/255);
reT=resultT/255;
end