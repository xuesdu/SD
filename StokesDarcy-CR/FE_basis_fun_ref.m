function out = FE_basis_fun_ref(x,y,basis_type,basis_index,basis_der_x,basis_der_y)

if basis_type == 200   % P0 element
    out = 1; 
elseif basis_type == 201 % CR 
    if basis_index == 1
        if basis_der_x == 0 && basis_der_y == 0
            out = 2*x+2*y-1;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 2;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 2;
        end
    elseif basis_index == 2
        if basis_der_x == 0 && basis_der_y == 0
            out = 1-2*x;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = -2;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 0;
        end
    elseif basis_index == 3
        if basis_der_x == 0 && basis_der_y == 0
            out = 1-2*y;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 0;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = -2;
        end
    end
end
        