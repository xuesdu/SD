function plot_cover( G, xcoord, ycoord, cover, edge_weight_flag)
% plot path covers of graph
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

if (edge_weight_flag == 0)
    pic = plot(G, 'XData', xcoord, 'YData', ycoord, 'LineWidth', 3);
else
    pic = plot(G, 'XData', xcoord, 'YData', ycoord, 'EdgeLabel',G.Edges.Weight);
end

num_cover = length(cover);

for i=1:num_cover
    highlight(pic,cover{i},'NodeColor','r','EdgeColor','r','LineWidth', 8);
end


end

