function out = FE_basis_fun_local(x,y,vertices,basis_type,basis_index,basis_der_x,basis_der_y)

if basis_type == 100
    out = 1;
elseif basis_type == 200  
    out = 1;
elseif basis_type == 201
    xn1 = vertices(1,1);  yn1 = vertices(2,1);
    xn2 = vertices(1,2);  yn2 = vertices(2,2);
    xn3 = vertices(1,3);  yn3 = vertices(2,3);
    J11 = xn2-xn1; J12 = xn3-xn1;
    J21 = yn2-yn1; J22 = yn3-yn1;
    J = J11*J22 - J12*J21;
    xh = ((yn3-yn1)*(x-xn1) - (xn3-xn1)*(y-yn1))/J;
    yh = (-(yn2-yn1)*(x-xn1) + (xn2-xn1)*(y-yn1))/J;
    
    if basis_der_x == 0 && basis_der_y == 0
        out =  FE_basis_fun_ref(xh,yh,basis_type,basis_index,basis_der_x,basis_der_y);
    elseif basis_der_x == 1 && basis_der_y == 0
        out = J22/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,1,0) +...
            -J21/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,0,1);
    elseif basis_der_x == 0 && basis_der_y == 1
        out = -J12/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,1,0)+...
            J11/J*FE_basis_fun_ref(xh,yh,basis_type,basis_index,0,1);
    end
end

