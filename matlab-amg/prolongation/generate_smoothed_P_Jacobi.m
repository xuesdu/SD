function [ P ] = generate_smoothed_P_Jacobi( tentative_P, As )
% Construct smoothed prolongation P using Jacobi method
%
% @ Xiaozhe Hu, Tufts University

%-------------------
% local variables
%-------------------
% size of A 
N = size(As,1);

smooth_weight = 0.67;
Dinv = spdiags(1./diag(As),0,N,N);

%-------------------
% smooth the tentative prolongation
%-------------------
%P = (speye(N) - smooth_weight*Dinv*As)*tentative_P;
P = tentative_P - smooth_weight*(Dinv*(As*tentative_P));

end