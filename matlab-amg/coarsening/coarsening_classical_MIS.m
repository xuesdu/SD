function [isC] = coarsening_classical_MIS(A, MIS_radius)
% Coarsening based on MIS and get CF-splitting for classical AMG
%
% @ Xiaozhe Hu, Tufts University 

%-----------------
% form power of A
%-----------------
A_radius = A;
for i = 1:MIS_radius-1
    A_radius = A_radius*A;
end

%-----------------
% find MIS 
%-----------------
%isC = MIS(A_radius);
isC = greedy_MIS(A_radius);

end

