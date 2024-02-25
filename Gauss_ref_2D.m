function [Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn)

if Gpn == 3
    Gauss_weights_ref_2D = [1/6 1/6 1/6];
    Gauss_nodes_ref_2D = [1/2 0;1/2 1/2;0 1/2];
elseif Gpn == 4
    Gauss_weights_ref_2D=[(1-1/sqrt(3))/8,(1-1/sqrt(3))/8,(1+1/sqrt(3))/8,(1+1/sqrt(3))/8];
    Gauss_nodes_ref_2D=[(1/sqrt(3)+1)/2,(1-1/sqrt(3))*(1+1/sqrt(3))/4;(1/sqrt(3)+1)/2,(1-1/sqrt(3))*(1-1/sqrt(3))/4;(-1/sqrt(3)+1)/2,(1+1/sqrt(3))*(1+1/sqrt(3))/4;(-1/sqrt(3)+1)/2,(1+1/sqrt(3))*(1-1/sqrt(3))/4];
elseif Gpn == 9
    Gauss_weights_ref_2D=[64/81*(1-0)/8,100/324*(1-sqrt(3/5))/8,100/324*(1-sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,40/81*(1-0)/8,40/81*(1-0)/8,40/81*(1-sqrt(3/5))/8,40/81*(1+sqrt(3/5))/8];
    Gauss_nodes_ref_2D=[(1+0)/2,(1-0)*(1+0)/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1+sqrt(3/5))/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1-sqrt(3/5))/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1+sqrt(3/5))/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1-sqrt(3/5))/4;(1+0)/2,(1-0)*(1+sqrt(3/5))/4;(1+0)/2,(1-0)*(1-sqrt(3/5))/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1+0)/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1+0)/4];
else
    return
end
    