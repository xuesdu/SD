function [ seed ] = Barabasi_network( Nodes,mlinks,num_seeds )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

initial_nodes = floor(Nodes/num_seeds);
seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0];

for i = 1:num_seeds
    T = tic;
    seed = SFNG(initial_nodes*i, mlinks, seed);
    T = toc(T);
    string = sprintf('network %d finished in %f seconds\n',i,T);
    disp(string);
end



end

