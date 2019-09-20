function b1 = binary_part(gray_image,a,b)
% ���þֲ���ֵ���������ֲ�����otus��
[m,n]=size(gray_image); 
result=zeros(m,n); 
for i=1:a:m 
    for j=1:b:n 
        if ((i+a)>m)&&((j+b)>n)     %�ֿ� 
            block1=gray_image(i:end,j:end); %���½�����
        elseif ((i+a)>m)&&((j+b)<=n) 
            block1=gray_image(i:end,j:j+b-1); %���������� 
        elseif ((i+a)<=m)&&((j+b)>n) 
            block1=gray_image(i:i+a-1,j:end);  %����������
        else 
             block1=gray_image(i:i+a-1,j:j+b-1); %��ͨ����
        end 
       
       [ block,~]=binary_otus(block1);
        if ((i+a)>m)&&((j+b)>n)            %�ϲ���� 
            result(i:end,j:end)=block; 
        elseif ((i+a)>m)&&((j+b)<=n) 
           result(i:end,j:j+b-1)=block; 
        elseif ((i+a)<=m)&&((j+b)>n) 
            result(i:i+a-1,j:end)=block; 
        else 
            result(i:i+a-1,j:j+b-1)=block; 
        end 
    end 
end 

b1=result;
end