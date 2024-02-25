function out = FE_basis_fun_local(x,y,vertices,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y)

xn1 = vertices(1,1);  yn1 = vertices(2,1);
xn2 = vertices(1,2);  yn2 = vertices(2,2);
xn3 = vertices(1,3);  yn3 = vertices(2,3);
J = (xn2-xn1)*(yn3-yn1) - (xn3-xn1)*(yn2-yn1);
xh = ((yn3-yn1)*(x-xn1) - (xn3-xn1)*(y-yn1))/J;
yh = (-(yn2-yn1)*(x-xn1) + (xn2-xn1)*(y-yn1))/J;

if basis_type == 100
    if basis_der_x == 0 && basis_der_y == 0
        out = 1;
    else
        out = 0;
    end
elseif basis_type == 200
    if basis_der_x == 0 && basis_der_y == 0
        out = 1;
    else
        out = 0;
    end
else
    if basis_der_x == 0 && basis_der_y == 0
        out = FE_basis_fun_ref(xh,yh,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y);
    elseif basis_der_x == 1 && basis_der_y == 0
        out = (yn3-yn1)/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,basis_vector,1,0) +...
            (yn1-yn2)/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,basis_vector,0,1);
    elseif basis_der_x == 0 && basis_der_y == 1
        out = (xn1-xn3)/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,basis_vector,1,0)+...
            (xn2-xn1)/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,basis_vector,0,1);
    end
end