function [ z ]=AMG_prod_prec(r, level, save_amgData, amgParam)
% AMG preconditioner (product of hierarchy)
%
% @ Xiaozhe Hu, Junyuan (Joanne) Lin, Tufts University, Ludmil Zikatanov,
% Penn State

%-------------------
% Preparation 
%-------------------
% max_it = amgParam.prec_max_it;

N = size(r,1);
z = zeros(N,1); 

%-------------------
% Main loop 
%-------------------
for k = 1:length(save_amgData)
    % call multigrid
    z = AMG_Cycle(save_amgData{k,1}, r, z, level, amgParam);
end

for k = length(save_amgData):-1:1
    % call multigrid
    z = AMG_Cycle(save_amgData{k,1}, r, z, level, amgParam);
end

end