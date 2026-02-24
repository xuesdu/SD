function [ x ] = cascade_eig(amgData, level, amgParam)
% cascadic Multigrid for Fiedler vector eigensolver
%
% X.Hu & J. Urschel

max_level = amgData(1).max_level;

n = amgData(level).N;

nsmooth = 2^(level-1)*amgParam.n_postsmooth;
smooth_type = amgParam.smooth_type;

if (level==max_level)
    
    m = amgData(1).number_eigen;
    I = speye(n, n);
    [V,~] = eigs(amgData(level).A + 0.1*I, m+1, 'SM');
    x = V(:,2:m+1);
    
else
    
    % coarse grid correction
    xH = cascade_eig(amgData, level+1, amgParam);
    
    % prolongation
    x = amgData(level).P*xH;
    
    % postsmoothing
    switch upper(smooth_type)
        case 'POWER'
            x = graphpowerit(x,amgData(level).A);
        case 'JACOBI'
            x = jacobi_GL(amgData(level).A, amgData(level).D, x, nsmooth);
        case 'FGS'
            x = forward_gs_GL(amgData(level).A, amgData(level).DL, x, nsmooth);
        case 'BGS'
            x = backward_gs_GL(amgData(level).A, amgData(level).DU, x, nsmooth);
        case 'SGS'
            x = symmetric_gs_GL(amgData(level).A, amgData(level).DL, amgData(level).DU, x, nsmooth);
    end
        
end

end % end of function 