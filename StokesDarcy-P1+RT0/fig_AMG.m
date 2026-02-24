% reponse to the reviewer

mu_vec = [1e-4, 1e-2, 1, 1e+2, 1e+4];
K_vec = [1e-4, 1e-2, 1, 1e+2, 1e+4];

%% use AMG directly
figure; hold on
% K = 1e-4
iterAMG = [14 17 12 12 13]; 
plot(mu_vec, iterAMG, 'ro-','MarkerFaceColor','r','LineWidth',2);
% K = 1e-2
iterAMG = [20 12 12 13 13]; 
plot(mu_vec, iterAMG, 'o-','color',[1, 0.5 0],'MarkerFaceColor',[1, 0.5 0],'LineWidth',2);
% K = 1
iterAMG = [12 12 13 13 13]; 
plot(mu_vec, iterAMG, 'o-', 'color',[0.3010 0.7450 0.9330], 'MarkerFacecolor', [0.3010 0.7450 0.9330],'LineWidth',2);
% K = 1e+2
iterAMG = [12 13 13 13 13]; 
plot(mu_vec, iterAMG, 'o-','color',[0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0],'LineWidth',2);
% K = 1e+4
iterAMG = [14 14 14 13 13]; 
plot(mu_vec, iterAMG, 'bo-','MarkerFaceColor','b','LineWidth',2);

set(gca, 'XScale', 'log');  
set(gca, 'XTick', mu_vec);
set(gca, 'XTickLabel', {'10^{-4}','10^{-2}','10^0', '10^2', '10^4'});
set(gca,'FontSize',14)
set(gcf,'Position',[300,300,500,450]); 
legend('$K = 10^{-4}$','$K = 10^{-2}$','$K = 10^0$', '$K = 10^2$', '$K = 10^4$','Interpreter','latex','FontSize',18);

xlabel('$\mu$', Interpreter='latex', FontSize=18);
ylabel('Number of iteration steps', 'Interpreter','latex','FontSize',18);
% title('Apply the UA-AMG directly')
box on;

%% outer iteration using GMRes-inner  
figure; hold on
% K = 1e-4
iterAMG = [9 9 9 10 10]; 
plot(mu_vec, iterAMG, 'ro-','MarkerFaceColor','r','LineWidth',2);
% K = 1e-2
iterAMG = [10 9 10 10 10]; 
plot(mu_vec, iterAMG, 'o-','color',[1, 0.5 0],'MarkerFaceColor',[1, 0.5 0],'LineWidth',2);
% K = 1
iterAMG = [8 8 10 10 10]; 
plot(mu_vec, iterAMG, 'o-', 'color',[0.3010 0.7450 0.9330], 'MarkerFacecolor', [0.3010 0.7450 0.9330],'LineWidth',2);
% K = 1e+2
iterAMG = [8 9 10 10 10]; 
plot(mu_vec, iterAMG, 'o-','color',[0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0],'LineWidth',2);
% K = 1e+4
iterAMG = [10 11 12 10 10]; 
plot(mu_vec, iterAMG, 'bo-','MarkerFaceColor','b','LineWidth',2);
ylim([4,20]);
set(gca, 'XScale', 'log');  
set(gca, 'XTick', mu_vec);
set(gca, 'XTickLabel', {'10^{-4}','10^{-2}','10^0', '10^2', '10^4'});
set(gca,'FontSize',14)
set(gcf,'Position',[300,300,500,450]); 
legend('$K = 10^{-4}$','$K = 10^{-2}$','$K = 10^0$', '$K = 10^2$', '$K = 10^4$','Interpreter','latex','FontSize',18);

xlabel('$\mu$', Interpreter='latex', FontSize=18);
ylabel('Number of iteration steps', 'Interpreter','latex','FontSize',18);
% title('Use GMRes as an inner solver')
box on;



%% The CPU time
figure; hold on
% K = 1e-4
iterAMG = [11.92 19.15 9.68 9.68 8.59]; 
plot(mu_vec, iterAMG, 'ro-','MarkerFaceColor','r','LineWidth',2);
% K = 1e-2
iterAMG = [26.51 9.44 9.98 8.45 7.74]; 
plot(mu_vec, iterAMG, 'o-','color',[1, 0.5 0],'MarkerFaceColor',[1, 0.5 0],'LineWidth',2);
% K = 1
iterAMG = [9.28 7.82 8.52 7.48 7.76]; 
plot(mu_vec, iterAMG, 'o-', 'color',[0.3010 0.7450 0.9330], 'MarkerFacecolor', [0.3010 0.7450 0.9330],'LineWidth',2);
% K = 1e+2
iterAMG = [12.55 6.90 7.87 7.57 7.62]; 
plot(mu_vec, iterAMG, 'o-','color',[0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0],'LineWidth',2);
% K = 1e+4
iterAMG = [11.53 7.53 7.73 7.79 7.91]; 
plot(mu_vec, iterAMG, 'bo-','MarkerFaceColor','b','LineWidth',2);

% cpu time for gmres
time = [7.52 13.13 7.62 7.80 6.13]; % gmres
plot(mu_vec, time, 'r+:','MarkerFaceColor','r','LineWidth',2);
time = [16.19 7.37 7.70 6.44 5.74];
plot(mu_vec, time, '+:','color',[1, 0.5 0],'MarkerFaceColor',[1, 0.5 0],'LineWidth',2);
time = [5.28 6.58 6.03 6.48 5.56];
plot(mu_vec, time, '+:', 'color',[0.3010 0.7450 0.9330], 'MarkerFacecolor', [0.3010 0.7450 0.9330],'LineWidth',2);
time = [4.32 5.86 6.14 5.77 6.06];
plot(mu_vec, time, '+:','color',[0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0],'LineWidth',2);
time = [6.53 9.83 10.92 5.86 5.93];
plot(mu_vec, time, '+:','MarkerFaceColor','b','LineWidth',2);


set(gca, 'XScale', 'log');  
set(gca, 'XTick', mu_vec);
set(gca, 'XTickLabel', {'10^{-4}','10^{-2}','10^0', '10^2', '10^4'});
set(gca,'FontSize',14)
set(gcf,'Position',[300,300,500,450]); 
legend('10^{-4}','10^{-2}','10^0', '10^2', '10^4','10^{-4}','10^{-2}','10^0', '10^2', '10^4');

xlabel('$\mu$', Interpreter='latex', FontSize=18);
ylabel('The CPU time (s)', 'Interpreter','latex','FontSize',18);
% title('Use GMRES as an inner solver')
box on;



