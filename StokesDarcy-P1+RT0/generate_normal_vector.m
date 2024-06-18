function normal = generate_normal_vector(vertices)

for k = 1:3
    if k == 1
        x1 = vertices(1,2); y1 = vertices(2,2);
        x2 = vertices(1,3); y2 = vertices(2,3);
        tau = [x1-x2;y1-y2];
        n = [tau(2);-tau(1)];
        if dot(n,[1;1])>0
            normal(:,1) = n;
        else
            normal(:,1) = -n;
        end
    elseif k == 2
        x1 = vertices(1,3); y1 = vertices(2,3);
        x2 = vertices(1,1); y2 = vertices(2,1);
        tau = [x1-x2;y1-y2];
        n = [tau(2);-tau(1)];
        if dot(n,[-1;0])>0
            normal(:,k) = n;
        else
            normal(:,k) = -n;
        end
    elseif k == 3
        x1 = vertices(1,1); y1 = vertices(2,1);
        x2 = vertices(1,2); y2 = vertices(2,2);
        tau = [x1-x2;y1-y2];
        n = [tau(2);-tau(1)];
        if dot(n,[0;-1])>0
            normal(:,k) = n;
        else
            normal(:,k) = -n;
        end
    end
end