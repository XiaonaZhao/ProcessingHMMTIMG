function b1 = binary_bernsen( f1 )
% �ֲ���ֵ���������ֲ�������ü���ֵ��
[m,n]=size(f1);
s=15;
t1=20;
exI=uint8(ones(m+2,n+2));%��չͼƬ��Ԥ����һ������
re=uint8(ones(m,n));
exI(2:m+1,2:n+1)=f1;%��ԭͼƬ��������
%==========�Ծ���������==========%
exI(1,:)=exI(2,:);
exI(m+2,:)=exI(m+1,:);
exI(:,1)=exI(:,2);
exI(:,n+2)=exI(:,n+1);

for i=2:m+1
    for j=2:n+1
        %===========��3*3�����ڵ���ֵ����ͼ����ж�ֵ�����������re��==========%
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