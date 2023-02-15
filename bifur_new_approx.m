%% Code for bifurcation plots of the local correlation approximation in the article "Finding analytical approximations for discrete, stochastic, individual-based models of ecology"

clear all
close all

%% r as a bifurcation parameter
T = 5000;


s=3;
T = 5000;

r_vals=1:1:30;

y1 = zeros(T, length(r_vals));

j=0;
for r = r_vals
    r
    j=j+1;
    y1(1,j) = rand(1);
  
  
    
    
    for t=1:T-1
        
        
        y1(t+1,j) = r*y1(t,j)*(1-1/(2*s+1)^2)^(r-1)*(y1(t,j)*(1-(1/(2*s+1)^2))^r-y1(t,j)+1)^((2*s+1)^2-1);
      
        
        
    end
    
    
end

figure;
plot(r_vals, (r_vals.*y1((4500:5000),:)),'.', 'color', [0, 0.4470, 0.7410], 'MarkerSize', 10)
% %plot(r_vals, y((2700:3001),:),'.', 'color', [0, 0.4470, 0.7410], 'MarkerSize', 15)
h=gca;
set(h,'box','off')
h.TickLength = [0.01 0.01];
h.LineWidth = 1.2;
 xticks([0 5 10 15 20 25 30])
 ylim([0 12])
 yticks([0 2 4 6 8 10 12])
set(h,'TickLabelInterpreter','latex')
ylabel('Population density', 'Interpreter', 'latex', 'FontSize', 12 )
xlabel('Number of offspring: $r$', 'Interpreter', 'latex', 'FontSize', 12)
set(h, 'FontSize', 14)
%title('()', 'Interpreter', 'latex')
set(get(gca,'title'),'Position',[15 11 0])
set(gca,'TickDir','out');
% % filename=sprintf('bifdiagram_s%d',d);
% % print(filename, '-depsc');
% % savefig(filename);


% %% s as bifrucation parameter
% s_vals = 1:100;
%    D=201;
% 
% 
% T = 5000;
% 
% r=30;
% 
% 
% y1 = zeros(T, length(s_vals));
% y2 = zeros(T, length(s_vals));
% 
% for s = s_vals
%     s
%     
%     y1(1,s) = rand(1);
%     y2(1,s) = rand(1);
%     
%     
%     delta = (s*(s+1))/(D*(2*s+1));
%     
%     
%     for t=1:T-1
%         
%         
%         y1(t+1,s) = (1-delta)*r*y1(t,s)*exp(-r*y1(t,s)) + delta*r*y2(t,s)*exp(-r*y2(t,s));
%         
%         y2(t+1,s) = (1-delta)*r*y2(t,s)*exp(-r*y2(t,s)) + delta*r*y1(t,s)*exp(-r*y1(t,s));
%         
%         
%     end
%     
%     
% end
% 
% 
% figure;
% plot(s_vals(1:14), 0.5*(r*y1((4500:5000),(1:14))+r*y2((4500:5000),(1:14))),'.', 'color', [0, 0.4470, 0.7410], 'MarkerSize', 5)
% hold on 
% plot(s_vals(71:100), 0.5*(r*y1((4500:5000),(71:100))+r*y2((4500:5000),(71:100))),'.', 'color', [0, 0.4470, 0.7410], 'MarkerSize', 5)
% hold on
% plot(s_vals(15:70), 0.5*(r*y1((4500:5000),(15:70))+r*y2((4500:5000),(15:70))),'.', 'color', [0, 0.4470, 0.7410], 'MarkerSize', 10)
% % %plot(r_vals, y((2700:3001),:),'.', 'color', [0, 0.4470, 0.7410], 'MarkerSize', 15)
% h=gca;
% set(h,'box','off')
% h.TickLength = [0.01 0.01];
% h.LineWidth = 1.2;
% xticks([0 20 40 60 80 100])
% % xticks([0 10 20 30])
%  ylim([0 12])
% yticks([0 2 4 6 8 10 12])
% set(h,'TickLabelInterpreter','latex')
% ylabel('Population density', 'Interpreter', 'latex', 'FontSize', 12 )
% xlabel('Dispersal distance: $s$', 'Interpreter', 'latex', 'FontSize', 12)
% set(h, 'FontSize', 24)
% %title('()', 'Interpreter', 'latex')
% set(get(gca,'title'),'Position',[15 11 0])
% set(h, 'FontSize', 14)
% % % filename=sprintf('bifdiagram_s%d',d);
% % % print(filename, '-depsc');
% % % savefig(filename);
