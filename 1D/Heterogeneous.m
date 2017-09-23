% CA for 1D with probablistic nucleation and heterogeneously deformed

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
        
        grain(1,i) = k;
        
    end 
       
end

mesh(grain);
view(2);
pause(0);

for i = 2:steps
    
    for j = 2:n-1
        
        left = cell(j-1,1);
        center = cell(j,1);
        right = cell(j+1,1);
        
        if (center == 1) 
            
            cell1(j,1) = 1;
            
            grain(i,j) = grain(i-1,j);        
        
        elseif (left == 0) && (right == 0)
            
            cell1(j,1) = 0;
            
            grain(i,j) = 0;
                       
        elseif (left == 1) 
            
            r2 = rand();
            
            if r2 < 0.4
                            
                cell1(j,1) = 1;        
            
                grain(i,j) = grain(i-1,j-1);
                
            else
                
                cell1(j,1) = 0;
                
                grain(i,j) = 0;
                
            end
            
        elseif (right == 1)
            
            r2 = rand();
            
            if r2 < 0.4
            
                cell1(j,1) = 1;
            
                grain(i,j) = grain(i-1,j+1); 
                
            else
                
                cell1(j,1) = 0;
                
                grain(i,j) = 0;
            
            end
      
        end
             
    end
    
    cell = cell1;
   
mesh(grain);
view(2);
pause(0); 
end
       
        
        
        
