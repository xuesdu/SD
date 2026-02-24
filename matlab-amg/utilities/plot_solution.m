function plot_solution( solu )
% plot solution on a N times N mesh
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

N = sqrt(length(solu));
mat_solu = reshape(solu, [N,N]);

mesh(1:N, 1:N, mat_solu);

% if (level_set_flag>0)
%     hold on;
%     contour(mat_solu);
% end


end

