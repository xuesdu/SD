function [ z ] = disagg_inexact_scale_prec(r, LS_dis, Ds, Ps, amgData, amgParam)
% Disaggregation preconditioner (inexact solve by AMG on the disaggregate graph)
%
% @ Xiaozhe Hu, Tufts University

% get size
N = size(LS_dis,1);

% AMG parameters
level = 1;
max_it = amgParam.prec_max_it;

% initialize
z_dis = zeros(N,1); 

% preconditioning
r_dis = Ps*r;  

% solve by AMG
for k = 1:max_it
    
    % call multigrid
    z_dis = AMG_Cycle(amgData, r_dis, z_dis, level, amgParam);

end

z =  Ps'*z_dis;

end
