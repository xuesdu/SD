function [ aggregation, num_agg ]=coarsening_agg_HEC_special( L )
% Coarsening based on heavy edge coarsening
% special version that some nodes are treated in a special way

n = length(L(:,1));
num_agg = 0;

%p=symrcm(L);
%p = amd(L);
%p = 1:n;
p = randperm(n);

aggregation=zeros(n,1);

% main loop
for i=1:n
    
    if aggregation(p(i))==0
    
        [~,m]=min(L(:,p(i)));
        
        if aggregation(m)==0
        
            num_agg = num_agg+1;
            aggregation(m) = num_agg;
            aggregation(p(i)) = num_agg;
        
        else
            
            aggregation(p(i)) = aggregation(m);
        
        end
        
    end
    
end

% keep the last 12 rigid body motion
num_agg = num_agg + 12;
aggregation(end-11:end) = (num_agg-11:num_agg);

     
        
