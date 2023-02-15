%% Code for extinction plots for the article "Finding analytical approximations for discrete, stochastic, individual-based models of ecology"


clear prob_rd

r_vals=0:0.05:40;
D_vals=0:0.05:201;
i=0;
for r=r_vals
    i=i+1;
    j=0;
for D= D_vals
   j=j+1;
prob_rd(i, j)= (1-r*exp(-1).*exp(-r*exp(-1))).^(D*D);
end 
end



figure; 
imagesc(D_vals, r_vals, prob_rd)
hold on 
%plot(D_vals, exp(1)*log(D_vals.^2), 'k')
%plot(D_vals, 2.8*log(D_vals.^2), 'k')
%plot(extinct2(:,1), extinct2(:,2),'xk', 'LineWidth',1.5)
hold off 

h = gca;
h.YDir = 'normal';


set(h,'box','on')
xlim([0 200])
ylim([0 40])
h.LineWidth =1.0;
set(h,'TickLabelInterpreter','latex')
ylabel('Number of offspring: $r$', 'Interpreter', 'latex', 'FontSize', 20 )
xlabel('Lattice size: $D$', 'Interpreter', 'latex', 'FontSize', 20)
set(h, 'FontSize', 12)
c = colorbar;
set(c, 'TickLabelInterpreter', 'latex');
set(c,'YTick',0:0.2:1);
set(gca,'TickDir','out');
h.TickLength = [0.01 0.0];
set(h,'box','off')
xlim([0 200])
c.LineWidth = 0.8;
hold on
errorbar(simulation_vals(:,1),simulation_vals(:,2), simulation_vals(:,3),"k",'LineWidth',0.6)
xticks([0 50 100 150 200])
% r = 1:1:40;
% d = 1:1:200;
% [D,R] = meshgrid(d,r);
% Z = (1-R.*exp(-1).*exp(-R.*exp(-1))).^(D.*D);
% % figure;
% % imagesc(Z)
% hold on
% contour(D,R,Z,'ShowText','off', 'color', 'k')
% %clabel(C,h)


% timesteps = 500; %number of timesteps
% 
% prob_less = prob_rd<1/timesteps;
% 
% figure;
% imagesc(D_vals, r_vals, prob_less)
% ax = gca;
% ax.YDir = 'normal';
% 
% lines = zeros(r,1);
% 
% 
% for j = 1:26
% 
%    hej(j) = find(prob_less(j,:), 1, 'first');
% end
