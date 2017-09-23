% CA for 1D with probablistic nucleation and single grain

clear all
clc

n = 1000;
steps = 200;
p_nucl = 0.01;

cell = zeros(n,1);
cell1 = zeros(n,1);

grain = zeros(steps,n);

r = rand(n,1);
r1 = n*r + 1;

n1 = n - n*p_nucl;

k = 0;

for i = 1:n
    
    if r1(i,1) >= n1
        
        k = k + 1;
        
        cell(i,1) = 1;
        
        grain(1,i) = 1;
        
    end 
       
end

mesh(grain);
view(2);
pause(0);

left = 0;

for i = 2:steps
    
    for j = 1:n-1
        
        center = cell(j,1);
        right = cell(j+1,1);
        
        if (left == 0) && (right == 0) && (center == 0)
            
            cell1(j,1) = 0;
            
        elseif (left == 1) || (right == 1) || (center == 1)
            
            cell1(j,1) = 1;        
            
        end
        
        left = cell(j,1);           
        
    end
        
    cell = cell1;
        
    if cell(n-1,1) == 1
        
        cell(n,1) = 1;
            
    end

    for j = 1:n
        
        grain(i,j) = cell(j,1);
        
    end
    
mesh(grain);
view(2);
pause(0); 
end
       
        
        
        
