function [ A ] = assembleConvectionDiffusion(N, epsilon)
% assemble upwind finite difference method for solving  
% - \varepsilon \Delta u + u_x = f
% u = 0 on boundary
% on unit square using uniform grid
% 
% Copyright (C)  Xiaozhe Hu.

% mesh size
h = 1/(N+1);

% some useful matrix
e = ones(N,1);
I = speye(N,N);

% assemble the diffusion part
D1d = spdiags([-1*e 2*e -1*e], -1:1, N, N);
Diffusion =  (kron(D1d, I) + kron(I, D1d));

% assemble the convection part
C1d = spdiags([-1*e 1*e], -1:0, N, N);
Convection = kron(I, C1d);
%Convection = kron(C1d, I);

% put two parts together
%A = (epsilon * Diffusion)/(h^2) + Convection/h;
A = (epsilon * Diffusion) + Convection*h;

end