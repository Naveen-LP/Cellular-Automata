%2D recrystallization with Moore boundary condtion

clear all
clc

nx = 500;
ny = 500;
steps = 100;

prob_nucl = 0.001;

cell = zeros(nx,ny);
cell1 = zeros(nx,ny);

grain = zeros(nx,ny,steps);
grain1 = zeros(nx,ny);

k = 0;
r = rand(nx,ny);
   
for i = 1:nx
    for j = 1:ny
                    
        if r(i,j) < prob_nucl
                          
            cell(i,j) = 1;
            
            k = k + 1;
            grain(i,j,1) = k;
                
        end
        
        grain1(i,j) = grain(i,j,1);
            
     end
end

mesh(grain1);
view(2);
pause(0);

for x = 2:steps
    
    for i = 2:nx-1
        for j = 2:ny-1
                 
              if cell(i,j) == 1
                    
                    cell1(i,j) = 1;
                    
                    grain(i,j,x) = grain(i,j,x-1);
                    
                else if cell(i-1,j) == 1
                        
                        cell1(i,j) = 1;
                        
                        grain(i,j,x) = grain(i-1,j,x-1);
                        
                    else if cell(i+1,j) == 1
                            
                            cell1(i,j) = 1;
                            
                            grain(i,j,x) = grain(i+1,j,x-1);
                            
                        else if cell(i,j-1) == 1
                                
                                cell1(i,j) = 1;
                                
                                grain(i,j,x) = grain(i,j-1,x-1);
                                
                            else if cell(i,j+1) == 1
                                    
                                    cell1(i,j) = 1;
                                    
                                    grain(i,j,x) = grain(i,j+1,x-1);
                                    
                                else if cell(i-1,j-1) == 1
                                        
                                        cell1(i,j) = 1;
                                        
                                        grain(i,j,x) = grain(i-1,j-1,x-1);
                                        
                                    else if cell(i+1,j-1) == 1
                                            
                                            cell1(i,j) = 1;
                                            
                                            grain(i,j,x) = grain(i+1,j-1,x-1);
                                            
                                        else if cell(i-1,j+1) == 1
                                                
                                                cell1(i,j) = 1;
                                                
                                                grain(i,j,x) = grain(i-1,j+1,x-1);
                                                
                                            else if cell(i+1,j+1) == 1
                                                    
                                                    cell1(i,j) = 1;
                                                    
                                                    grain(i,j,x) = grain(i+1,j+1,x-1);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            
            grain1(i,j) = grain(i,j,x);
            
        end
     end
    
    cell = cell1;
    
    mesh(grain1);
    view(2);
    pause(0);
    
end

