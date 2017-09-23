% CA for 1D with nucleation in the center

clear all
clc

n = 1000;
steps = 200;

cell = zeros(n,1);
cell1 = zeros(n,1);

grain = zeros(steps,n);

cell(500,1) = 1;

left = 0;

for i = 1:steps
    
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
