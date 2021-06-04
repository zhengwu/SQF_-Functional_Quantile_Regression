%code for 1D FA curve analysis

clear all;
close all;


%% load data
load('data/ADNI_1D.mat')

%% orgnize X and Y;

%  Y
Y = CC_FA_Selected';
m = size(Y,2);
n = size(Y,1);

%  X
temp = COV_Selected';
clear X;
for i=1:n
   X(:,i) = [1;temp(:,i)];
end
X(2,:) = X(2,:) -1;
for i=3:size(X,1)
   X(i,:) = (X(i,:)  - mean(X(i,:)))/std(X(i,:));
end
X(4,:) = X(5,:);
X(5,:) = [];

%parameters for model fitting 
s = (0:m-1)/(m-1);
p = size(X,1);


h = 0.1;
a = 1;
Sigma = zeros(m,m);
for i=1:m
    for j=1:m
        
        Sigma(i,j) = a*exp(-((s(i)-s(j))/h)^2 );
        
    end
end



%% lambda selection - do not skip for new data
% tic
% tau = 0.5;
% lama  = [1:2:60,80:10:160,170:20:300,500:20:600,1000,2000];
% for i=1:length(lama)
%     lam=lama(i);
%     [bta, dd_m, Yhat, gacv(i)] = quan_PrimDual2(Y, X, Sigma, tau, lam);
%     
%     %plot true beta and fitted beta
%     figure(10),clf;
%     plot(s,bta,'--b');
%
%     pause(0.1);
% end
% toc
%find the best lambda
%lam = lama(find(gacv==min(gacv)));



taua = 0.01:0.01:0.99;
lam=300;
idx = 1;
for i=1:length(taua)
    tau = taua(i);
    [bta, dd_m, Yhat, gacv(i)] = quan_PrimDual2(Y, X, Sigma, tau, lam);
    betall(i,:,:) = bta;
end


%plot betas

colora = parula(99);

x_axis = [0 1 -0.2 0.2];
figure(10);clf;
for i=1:length(taua)
    if(mod(i,1)==0)
        bet_1 = squeeze(betall(i,1,:));
        figure(10);
        hold on;
        if(i==1)
            plot(s,bet_1,'--b','linewidth',3);
        elseif(i==99)
            plot(s,bet_1,'-.g','linewidth',3);
        elseif(i==50)
            plot(s,bet_1,'r','linewidth',3);
        else
            plot(s,bet_1,'linewidth',0.5,'color',colora(i,:));
        end
    end
end
set(gca,'fontsize',22);
%legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9');
title('Intercept')


figure(11);clf;
figure(11);clf;
for i=1:length(taua)
    if(mod(i,1)==0)
        bet_1 = squeeze(betall(i,2,:));
        figure(11);
        hold on;
        if(i==1)
            plot(s,bet_1,'--b','linewidth',3);
        elseif(i==99)
            plot(s,bet_1,'-.g','linewidth',3);
        elseif(i==50)
            plot(s,bet_1,'r','linewidth',3);
        else
            plot(s,bet_1,'linewidth',0.5,'color',colora(i,:));
        end
    end
end
hold on, plot(s,0*ones(length(s),1),'k','linewidth',1);
set(gca,'fontsize',22);
axis(x_axis)
%legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9');
title('Gender')

figure(12);clf;
for i=1:length(taua)
    if(mod(i,1)==0)
        bet_1 = squeeze(betall(i,3,:));
        figure(12);
        hold on;
        if(i==1)
            plot(s,bet_1,'--b','linewidth',3);
        elseif(i==99)
            plot(s,bet_1,'-.g','linewidth',3);
        elseif(i==50)
            plot(s,bet_1,'r','linewidth',3);
        else
            plot(s,bet_1,'linewidth',0.5,'color',colora(i,:));
        end
    end
end
hold on, plot(s,0*ones(length(s),1),'k','linewidth',1);
set(gca,'fontsize',22);
axis(x_axis)
%legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9');
title('Age')

figure(13);clf;
for i=1:length(taua)
    if(mod(i,1)==0)
        bet_1 = squeeze(betall(i,4,:));
        figure(13);
        hold on;
        if(i==1)
            plot(s,bet_1,'--b','linewidth',3);
        elseif(i==99)
            plot(s,bet_1,'-.g','linewidth',3);
        elseif(i==50)
            plot(s,bet_1,'r','linewidth',3);
        else
            plot(s,bet_1,'linewidth',0.5,'color',colora(i,:));
        end
    end
end
hold on, plot(s,0*ones(length(s),1),'k','linewidth',1);
set(gca,'fontsize',22);
axis(x_axis)
%legend([1,50,90],'0.1','0.5','0.9');
title('ADAS')