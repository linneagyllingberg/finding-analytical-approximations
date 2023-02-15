%% CODE FOR INDIVIDUAL BASED MODEL IN THE PAPER "Finding analytical approximations for discrete, stochastic, individual-based models of ecology"

% authors: Linnéa Gyllingberg, Åke Brönnström and David Sumpter

clear all
close all 


num_vals = 1:10;
for num_sim = num_vals
   
s_vals=51:100;


for s=s_vals
    s
%Size of the world 
D=2*s+1;
%Number of offspring

%Increase in fecundity from help.
%Maximum set to 240
help=0;,
maxhelp=240;

makemovie=0;
clear pop_take
T=5000;
pop_take = zeros(35,T);

 for b=1:19
     pop_take(b,:) = ones(1, 5000);
     
 end

for b=20:35
%Number of time steps

b;

disp([num_sim,s,b])
%Maximum dispersal distance
%s=1;

%Initial population
den=0.3;

%Stores the ID of the individual
a=zeros(T,D,D);
%Stores the strategy of individual
%-1 indicates that site has been occupied first by a type 2 (selfish)
%and no-one can use it.
c=zeros(T,D,D);
r=zeros(T,D,D);
%Baseline reproduction
effb=ones(T,D,D)*b;

counta=0;
for x=1:D
    for y=1:D
        if rand<den/2
            counta=counta+1;
            a(1,x,y)=counta;
            c(1,x,y)=2;
        end
        if rand<den/2
            counta=counta+1;
            a(1,x,y)=counta;
            c(1,x,y)=2;
        end
    end
end

%shiftdim(a(1,:,:));

for t=1:T-1
    
    %I randomly permute the order individuals are moved
    %This ensures that the first individual arriving at a site is randomly
    %chosen.

    %Stores the number arriving at a site, but not the first to arrive.
    
    numbercomingto=zeros(2,D,D);
    %Counts the number which are identical by decesent.
    related=zeros(D,D);
    
    
    
    for nextsite=randperm(D*D)
        x=rem(nextsite,D);
        if (x==0)
            x=D;
        end
        y=ceil(nextsite/D);
            if c(t,x,y)>0
                %b is at most twice its baseline rate.
                bhere=min(effb(t,x,y),maxhelp);
                
                %Site is occupied 
                smovex=floor(rand(bhere,1)*(2*s+1))-s;
                smovey=floor(rand(bhere,1)*(2*s+1))-s;
                newx=x+smovex;
                newy=y+smovey;

                %Calculate which have gone over the edge and give them all new Ids
                newNumber=((newx<1)+(newx>D)+(newy<1)+(newy>D))>0;
                
                
                
                if sum(newNumber)>0
                    %Give each in over the edge individual own ID
                    %counta=counta+b;
                    %moveIds=a(t,x,y).*(1-newNumber)+newNumber.*[counta+1:counta+b]';
%                   %Give each coming in over the edge individual same ID
                    counta=counta+1;
                    moveIds=a(t,x,y).*(1-newNumber)+newNumber*counta;
                else
                    moveIds=a(t,x,y)*ones(bhere,1);
                end
                moveStrategy=c(t,x,y);
                
                
                %Now make sure everything lands in the correct boxes
                newx=mod(newx,D)+D*(mod(newx,D)==0);
                newy=mod(newy,D)+D*(mod(newy,D)==0);

                %Move them to a new site.
                for i=1:bhere
                    if (c(t+1,newx(i),newy(i))~=0)
                        %Site is occupied, increment number going to that site,
                        %but the original one there owns the site (unless removed).
                        
                        %Calculate relatedness. First increment the number
                        %arriving at the site.
                        numbercomingto(moveStrategy,newx(i),newy(i))=numbercomingto(moveStrategy,newx(i),newy(i))+1;
                        %See if they are same by descent
                        if moveIds(i)==a(t+1,newx(i),newy(i))
                                related(newx(i),newy(i))=related(newx(i),newy(i))+1;
                        end
                        if (moveStrategy==1)
                            %Help out by increasing number of offspring.
                            effb(t+1,newx(i),newy(i))=effb(t+1,newx(i),newy(i))+help;                            
                        elseif (moveStrategy==2)
                            %Replace occupant if occupant is 1. Otherwise
                            %die
                            if ((c(t+1,newx(i),newy(i))==-1) || (c(t+1,newx(i),newy(i))==2))
                                c(t+1,newx(i),newy(i))=(-1);
                                %c=-1 means that a stategy 2 has already been at
                                %this site.
                                %Leave id of first individual there so we
                                %can measure relatednesss
                            else
                                a(t+1,newx(i),newy(i))=moveIds(i);
                                c(t+1,newx(i),newy(i))=2;
                                %Reset relatedness because all individuals
                                %already there are strategy 1 and thus not identical by
                                %desent.
                                related(newx(i),newy(i))=0;
                            end
                        end
                    else
                        %This one is there first it owns it and reproduces
                        %(unless replaced)
                        a(t+1,newx(i),newy(i))=moveIds(i);
                        c(t+1,newx(i),newy(i))=moveStrategy;
                    end
                end
            end
        
    end
    %shiftdim(a(t+1,:,:))
    
    %Relatedness. You ad yourself to this point
    r(t+1,:,:)=(related)./(max(1,shiftdim(sum(numbercomingto,1)))).*shiftdim(c(t+1,:,:)~=0);
    %Average relatedness.
    mr(b,t+1)=sum(sum(r(t+1,:,:)))./sum(sum(c(t+1,:,:)~=0));
    
    %For each independently
    r1=(related)./(max(1,shiftdim(sum(numbercomingto,1)))).*shiftdim(c(t+1,:,:)==1);
    mr1(b,t+1)=sum(sum(r1))./sum(sum(c(t+1,:,:)==1));
    r2=(related)./(max(1,shiftdim(sum(numbercomingto,1)))).*(shiftdim(c(t+1,:,:)==2)+shiftdim(c(t+1,:,:)==-1));
    mr2(b,t+1)=sum(sum(r2))./sum(sum((shiftdim(c(t+1,:,:)==2)+shiftdim(c(t+1,:,:)==-1))));
    
    
    
    if makemovie
        figure(1)
        imagesc(max(shiftdim(c(t,:,:)),0),[0 2])
        M(t)=getframe;
%         figure(2)
%         imagesc(max(shiftdim(r(t,:,:)),0),[0 1])
%         R(t)=getframe;
    end
end

%pop_help(b,:)=sum(sum(c(0:T,:,:)==1,2),3)';
 pop_take(b,:)=sum(sum(c(1:T,:,:)==2,2),3)';
% pop_take(b,:)=sum(sum(c(round(T/2):T,:,:)==2,2),3)';
% filename=sprintf('s_1_outputnonhelpers_new%d',b);
% save(filename)

end

% figure(8)
% subplot(3,3,1)
%     imagesc(max(shiftdim(c(t+1,:,:)),0)>0,[0 2])
% subplot(3,3,2)
%     imagesc(max(shiftdim(r(t+1,:,:)),0)>0,[0 1])
% subplot(3,3,3)    
%     imagesc(max(shiftdim(a(t+1,:,:)),0)>0,[0 2])
% subplot(3,3,4)
% imagesc((max(shiftdim(a(t+1,:,:)),0)>0)+(max(shiftdim(c(t+1,:,:)),0)>0),[0 2])
% subplot(3,3,5)
% imagesc(related,[0 1])
% subplot(3,3,6)
% imagesc(shiftdim(sum(numbercomingto,1)),[0 1])
% 
% 
% figure(9)
% for b=10:10:120
% subplot(4,3,b/10)
% hold off
% plot(mr(b,:),'b')
% hold on
% plot(mr1(b,:),'r')
% plot(mr2(b,:),'g')
% end
% 
% 
% 
% figure(11)
% hold off
% plot([1:120],mean(mr(:,round(T/2):T),2),'b')
% hold on
% plot([1:120],mean(mr1(:,round(T/2):T),2),'g')
% plot([1:120],mean(mr2(:,round(T/2):T),2),'r')
% hx=xlabel('Resource Capacity: b')
% set(hx,'FontSize',14)
% hy=ylabel('Relatedness (Probability of Identity by Descent')
% set(hy,'FontSize',14)
% 
% 
% figure(13)
% hb=mean(pop_help(:,round(T/2)),2)./(mean(pop_help(:,round(T/2)),2)+mean(pop_take(:,round(T/2)),2));
% rb=mean(mr(:,round(T/2):T),2);
% hlx=plot(rb,hb,'k-')
% hx=ylabel('Proportion of Helping Individuals')
% set(hx,'FontSize',14)
% hy=xlabel('Relatedness: Probability of Identity by Descent')
% set(hy,'FontSize',14)
% axis([-0.02 1.02 -0.02 1.02])
% 
% 
% figure(12)
% hold off
% plot([1:120],mean(pop_help(:,round(T/2)),2),'b')
% hold on
% plot([1:120],mean(pop_take(:,round(T/2)),2),'r')
% hx=xlabel('Resource Capacity: b')
% set(hx,'FontSize',14)
% hy=ylabel('Number of Individuals')
% set(hy,'FontSize',14)
% 
% aviobj = avifile('example_A100.avi');
% 
% for   t=1:1000
%         figure(1)
%         imagesc(max(shiftdim(c(t,:,:)),0),[0 2])
%         hx=xlabel('X')
%         set(hx,'FontSize',14)
%         hy=ylabel('Y')
%         set(hy,'FontSize',14)
%         M=getframe;
%         aviobj = addframe(aviobj,M);
% 
%         figure(2)
%         imagesc(max(shiftdim(r(t,:,:)),0),[0 1])
%         hx=xlabel('X')
%         set(hx,'FontSize',14)
%         hy=ylabel('Y')
%         set(hy,'FontSize',14)
%         R(t)=getframe;
% end
% aviobj = close(aviobj);
% 
%             
% for b=1:120
%     filename=sprintf('outputForbHelp%d',b);
%     load(filename)
%     pop_help(b,:)=sum(sum(c(round(T/2):T,:,:)==1,2),3)';
%     pop_take(b,:)=sum(sum(c(round(T/2):T,:,:)==2,2),3)';
%     save(filename)
% end
% 
% 
% 
%       
% 
% 

%%% Including bifurcation plot
b_vals=1:35;
for b=b_vals
    pop_take(b,:)= pop_take(b,:)*b;
end

filename2 = sprintf('new_D_equal_2splus1_poptake_s_%d_simnum_%d',s,num_sim);
save(filename2)

end


% for t=4900:5000
%     figure(1)
%     hold on
%     plot(b_vals,pop_take(:,t)./(D*D), 'k.')
%     xlabel('$r$', 'Interpreter','latex', 'FontSize', 14)
%     ylabel('Population density', 'Interpreter','latex', 'FontSize', 14)
% 
% end
% 
% filename=sprintf('Ds_same_IBM_bifur_s_%d',s);
% print(filename,'-depsc')
% savefig(filename)
% close all
end
% 
% 
% b_vals=1:8;
% for t=300:501
%     figure(3)
%     hold on
%     plot(b_vals,pop_help(:,t)./(D*D), 'k.')
%     xlabel('Reproductive rate, $r$', 'Interpreter','latex', 'FontSize', 14)
%     ylabel('Population density, $x$', 'Interpreter','latex', 'FontSize', 14)
% 
% end
% % print('IBM_s5_rto30','-depsc')
% % savefig('IBM_s5_rto30')