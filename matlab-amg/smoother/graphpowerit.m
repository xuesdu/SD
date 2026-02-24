% Graph Power Iteration
% X.Hu & J. Urschel
function [u2, k]=graphpowerit(u2,L)

%n=length(u2);
n = size(u2,1);
m = size(u2,2);
v1=ones(n,1)/sqrt(n);

%Defining Bg
g=max(sum(abs(L)));
Bg=g*speye(n,n)-L;

% get rid of constant
for j=1:m
    u2(:,j)=u2(:,j)-(u2(:,j)'*v1)*v1;
    u2(:,j)=u2(:,j)/norm(u2(:,j));
end
v2=zeros(n,1);

k = 0;

while v2'*u2 < 1-10^(-5)
    k = k+1;
    
    for j=1:m
        v2=u2(:,j);
        if sum(v2)>.01
            v2=v2-(v2'*v1)*v1;
        end
        u2(:,j)=Bg*v2;
        u2(:,j)=u2(:,j)/norm(u2(:,j));
    end
    
end

% find orthonormal basis
[u2, ~] = qr(u2, 0);

    



        
