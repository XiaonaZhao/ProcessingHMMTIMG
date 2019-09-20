function b1 = binary_bernsen( f1 )
% 局部二值化方法，局部区域采用简单阈值。
[m,n]=size(f1);
s=15;
t1=20;
exI=uint8(ones(m+2,n+2));%扩展图片，预分配一个矩阵
re=uint8(ones(m,n));
exI(2:m+1,2:n+1)=f1;%把原图片赋给矩阵
%==========对矩阵进行填充==========%
exI(1,:)=exI(2,:);
exI(m+2,:)=exI(m+1,:);
exI(:,1)=exI(:,2);
exI(:,n+2)=exI(:,n+1);

for i=2:m+1
    for j=2:n+1
        %===========求3*3区域内的阈值并对图像进行二值化，结果存在re中==========%
        ma=max(max(exI(i-1:i+1,j-1:j+1)));
        mi=min(min(exI(i-1:i+1,j-1:j+1)));
        t=(ma+mi)/2;
        if ma-mi>s
             if exI(i,j)>t
                 re(i-1,j-1)=255;
             else
                re(i-1,j-1)=0;
            end
        else 
            if t>t1
                re(i-1,j-1)=255;
             else
                re(i-1,j-1)=0;
            end
        end
        
    end
end

b1=re;
end