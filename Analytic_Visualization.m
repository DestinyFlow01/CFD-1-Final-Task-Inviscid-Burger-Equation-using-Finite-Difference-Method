%Visualizing the analytical solution of burger equation

clc; clear all; 
Nx = 1001; Ny = Nx;
u = zeros(Nx, Ny);
X = 1; Y = 1;
dx = X/(Nx-1); dy = Y/(Ny-1); 

%Coordinate matrix 
X = zeros(Nx, Ny);
Y = zeros(Nx, Ny);
for i = 1:Nx
    X(:,i) = (i-1)*dx;
end

for j = 1:Ny
    Y(Ny-j+1, :) = (j-1)*dy;
end

%Velocity solution : 
for i = 1:Nx
    for j = 1:Ny
        if(Y(i,j) <= 0.5)
            
            if(X(i,j) <= 1.5*Y(i,j))
                u(i,j) = 1.5;
            elseif(X(i,j) >= 1-0.5*Y(i,j))
                u(i,j) = -0.5;
            else 
                u(i,j) = (1.5-2*X(i,j))/(1-2*Y(i,j));
            end
                 
        elseif (Y(i,j)>=0.5)
            
            if(X(i,j) <= 0.5+0.5*Y(i,j))
                u(i,j) = 1.5;
            else
                u(i,j) = -0.5;
            end
                
        end
        
    end
end

%Visualization
a = surf(X,Y,u);
a.EdgeColor = 'none';
xlabel('X')
ylabel('Y')
colorbar