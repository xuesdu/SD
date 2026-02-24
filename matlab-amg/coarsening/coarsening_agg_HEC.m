function [ aggregation, num_agg ]=coarsening_agg_HEC( L )
% Coarsening based on heavy edge coarsening

n = length(L(:,1));
num_agg = 0;

%p=symrcm(L);
p = amd(L);
%p = 1:n;
%p = randperm(n);

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

     
        
