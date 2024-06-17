function normal = generate_normal_vector(edge_index,end_point_1,end_point_2)


tau1 = end_point_2(1,1) - end_point_1(1,1);
tau2 = end_point_2(2,1) - end_point_1(2,1);
n = [-tau2;tau1];
n_length = sqrt(tau1^2+tau2^2);
if edge_index == 1
    if dot(n,[1;1]) > 0
        normal = n;
    else
        normal = -n;
    end
elseif edge_index == 2
    if dot(n,[-1;0]) >0
        normal = n;
    else
        normal = -n;
    end
elseif edge_index == 3
    if dot(n,[0;-1])>0
        normal = n;
    else
        normal = -n;
    end
end
normal = normal/n_length;  % 单位化