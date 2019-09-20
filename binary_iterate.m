function b1 = binary_iterate(f1)
% 用迭代的方法对图像进行二值化

t = mean(f1(:));                  
is_done = false;                        
count = 0;                              
while ~is_done                         
    
    r1 = f1(f1 <= t);              
    r2 = f1(f1 > t); 
    temp1 = mean(r1(:));              
    if isnan(temp1)                  
         temp1 = 0; 
    end 
    temp2 = mean(r2(:));               
    if isnan(temp2)                      
         temp2 = 0; 
    end 
    t_new = (temp1 + temp2)/2; 
    is_done = abs(t_new - t) < 1;             
    t = t_new;                           
    count = count+1;                     
    if count >= 1000                    
        Error = 'Error:Cannot find the ideal threshold.';
        Return                       
    end 
end 
b1 = im2bw(mat2gray(f1), t/256);
end