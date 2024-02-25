function [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D)

xn1 = vertices(1,1); yn1 = vertices(2,1);  % 当前单元1号点的坐标(xn1,yn1)
xn2 = vertices(1,2); yn2 = vertices(2,2);
xn3 = vertices(1,3); yn3 = vertices(2,3);

J = abs((xn2-xn1)*(yn3-yn1) - (xn3-xn1)*(yn2-yn1));
Gauss_weights_2D = J*Gauss_weights_ref_2D;
Gauss_nodes_2D(:,1) = xn1+(xn2-xn1)*Gauss_nodes_ref_2D(:,1) + (xn3-xn1)*Gauss_nodes_ref_2D(:,2);
Gauss_nodes_2D(:,2) = yn1+(yn2-yn1)*Gauss_nodes_ref_2D(:,1) + (yn3-yn1)*Gauss_nodes_ref_2D(:,2);