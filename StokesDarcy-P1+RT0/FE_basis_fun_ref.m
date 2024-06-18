function out = FE_basis_fun_ref(x,y,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y)


if basis_type == 200 
    if basis_der_x == 0 && basis_der_y == 0
        out = 1;
    else
        out = 0;
    end
elseif basis_type == 201
    if basis_index == 1
        if basis_der_x == 0 && basis_der_y == 0
            out = 1-x-y;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = -1;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = -1;
        end
    elseif basis_index == 2
        if basis_der_x == 0 && basis_der_y == 0
            out = x;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 1;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 0;
        end
    elseif basis_index == 3
        if basis_der_x == 0 && basis_der_y == 0
            out = y;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 0;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 1;
        end
    end
elseif basis_type == 202
    if basis_index == 1
        if basis_der_x == 0 && basis_der_y == 0
            out = 2.*x.^2+2.*y.^2+4.*x.*y-3.*y-3.*x+1;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 4.*x+4.*y-3;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 4.*y+4.*x-3;
        elseif basis_der_x == 2 && basis_der_y == 0
            out = 4;
        elseif basis_der_x == 0 && basis_der_y == 2
            out = 4;
        elseif basis_der_x == 1 && basis_der_y == 1
            out = 4;
        end
    elseif basis_index == 2
        if basis_der_x == 0 && basis_der_y == 0
            out = 2.*x.^2-x;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 4.*x-1;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 0;
        elseif basis_der_x == 2 && basis_der_y == 0
            out = 4;
        elseif basis_der_x == 0 && basis_der_y == 2
            out = 0;
        elseif basis_der_x == 1 && basis_der_y == 1
            out = 0;
        end
    elseif basis_index == 3
        if basis_der_x == 0 && basis_der_y == 0
            out = 2.*y.^2-y;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 0;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 4.*y-1;
        elseif basis_der_x == 2 && basis_der_y == 0
            out = 0;
        elseif basis_der_x == 0 && basis_der_y == 2
            out = 4;
        elseif basis_der_x == 1 && basis_der_y == 1
            out = 0;
        end
    elseif basis_index == 4
        if basis_der_x == 0 && basis_der_y == 0
            out = -4.*x.^2-4.*x.*y+4.*x;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = -8.*x-4.*y+4;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = -4*x;
        elseif basis_der_x == 2 && basis_der_y == 0
            out = -8;
        elseif basis_der_x == 0 && basis_der_y == 2
            out = 0;
        elseif basis_der_x == 1 && basis_der_y == 1
            out = 0;
        end
    elseif basis_index == 5
        if basis_der_x == 0 && basis_der_y == 0
            out = 4.*x.*y;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = 4.*y;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = 4.*x;
        elseif basis_der_x == 2 && basis_der_y == 0
            out = 0;
        elseif basis_der_x == 0 && basis_der_y == 2
            out = 0;
        elseif basis_der_x == 1 && basis_der_y == 1
            out = 4;
        end
    elseif basis_index == 6
        if basis_der_x == 0 && basis_der_y == 0
            out = -4.*y.^2-4.*x.*y+4.*y;
        elseif basis_der_x == 1 && basis_der_y == 0
            out = -4.*y;
        elseif basis_der_x == 0 && basis_der_y == 1
            out = -8.*y-4.*x+4;
        elseif basis_der_x == 2 && basis_der_y == 0
            out = 0;
        elseif basis_der_x == 0 && basis_der_y == 2
            out = -8;
        elseif basis_der_x == 1 && basis_der_y == 1
            out = -4 ;
        end
    end
end


if basis_type == 2011
    if basis_index == 1
        if basis_vector == 1
            if basis_der_x == 0 && basis_der_y == 0
                out = sqrt(2)*x;
            elseif basis_der_x == 1 && basis_der_y == 0
                out = sqrt(2);
            elseif basis_der_x == 0 && basis_der_y == 1
                out = 0;
            end
        elseif basis_vector == 2
            if basis_der_x == 0 && basis_der_y == 0
                out = sqrt(2)*y;
            elseif basis_der_x == 1 && basis_der_y == 0
                out = 0;
            elseif basis_der_x == 0 && basis_der_y == 1
                out = sqrt(2);
            end
        end
    elseif basis_index == 2
        if basis_vector == 1
            if basis_der_x == 0 && basis_der_y == 0
                out = x-1;
            elseif basis_der_x == 1 && basis_der_y == 0
                out = 1;
            elseif basis_der_x == 0 && basis_der_y == 1
                out = 0;
            end
        elseif basis_vector == 2
            if basis_der_x == 0 && basis_der_y == 0
                out = y;
            elseif basis_der_x == 1 && basis_der_y == 0
                out = 0;
            elseif basis_der_x == 0 && basis_der_y == 1
                out = 1;
            end
        end
        
    elseif basis_index == 3
        if basis_vector == 1
            if basis_der_x == 0 && basis_der_y == 0
                out = x;
            elseif basis_der_x == 1 && basis_der_y == 0
                out = 1;
            elseif basis_der_x == 0 && basis_der_y == 1
                out = 0;
            end
        elseif basis_vector == 2
            if basis_der_x == 0 && basis_der_y == 0
                out = y-1;
            elseif basis_der_x == 1 && basis_der_y == 0
                out = 0;
            elseif basis_der_x == 0 && basis_der_y == 1
                out = 1;
            end
        end
    end
else
    return;
end