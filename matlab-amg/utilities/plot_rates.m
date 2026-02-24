function [] = plot_rates(nsize,cputimes1,cputimes2,cputimes3,cputimes4)
% make the plot of matrix size vs. cpu times with reference line
% cputime1-4 is cpu times for each method
% nsize is a vector of testing matrix sizes

% Junyuan Lin @ Tufts University
n = length(nsize);
ref_line = zeros(n,1);
ref_line2 = zeros(n,1);
ref_line(1) = (cputimes1(1))-2.3; % regular case's starting point
%ref_line2(1) = (cputimes1(1))-1.5;
%if size(cputimes,1)==2 % b = 0 case
    for i = 2:n
        ref_line(i) = ref_line(i-1)*4*log(4);
    end
    for i = 2:n
        ref_line2(i) = ref_line2(i-1)*8;
    end
%     for i = 2:n
%         rate1 = cputimes(1,i)/cputimes(1,i-1);
%     end
%     for i = 2:n
%         rate2 = cputimes(2,i)/cputimes(2,i-1);
%     end
    %plot(nsize,cputimes1,'-.r*', 'LineWidth',4)
    loglog(nsize,cputimes1,'-.r*', 'LineWidth',4)
    hold on
    %plot(nsize,cputimes2,'--bO','LineWidth',4)
    loglog(nsize,cputimes2,'-.bO','LineWidth',4)
    %hold on
    loglog(nsize,cputimes3,'-.g*', 'LineWidth',4)
    loglog(nsize,cputimes4,'-.m*', 'LineWidth',4)
    %plot(nsize,4*log(4)*ones(n,1))
    %plot(nsize,ref_line, 'k' ,'LineWidth',4)
    loglog(nsize,ref_line, '--k' ,'LineWidth',3)
    %loglog(nsize,ref_line2, ':r' ,'LineWidth',3)
    legend('PC-\alphaAMG(1) on low-f b','PC-\alphaAMG(2) on low-f b','PC-\alphaAMG(1) on zero-sum b','PC-\alphaAMG(2) on zero-sum b','slope = 1', 'Location','northwest')
    %legend('Algorithm 6','slope = 1', 'Location','northwest')
    xlabel('Matrix Size')
    ylabel('CPU time (s)')
    xlim([10^4, 1.5*10^6]);
%end
end

