function [Gauss_weights,Gauss_nodes] = Gauss_local_1D(lb,ub,Gauss_weights_ref_1D,Gauss_nodes_ref_1D)
% 通用code：lb: 左端点的坐标；ub: 右端点的坐标

Gauss_weights = (ub-lb)*Gauss_weights_ref_1D/2;
Gauss_nodes = (ub-lb)*Gauss_nodes_ref_1D/2 + (ub+lb)/2;