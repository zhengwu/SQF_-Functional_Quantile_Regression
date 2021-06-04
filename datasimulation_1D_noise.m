% data simulation for 1D FRQR
%by Zhengwu Zhang
%08/10/2019

% model:
% Q_Y(s_j)(tau) = mu_x(s_j) + sigma_x(s_j)Phi^{-1}(tau)

%clear all;
%close all;

dispidx = 0;
rand_indicator = 0;

% global parameters - for dimentionality of the functions
m = 100;
s = (0:m-1)/(m-1);
n = 500; % number of total function


%simulate X
for idx = 1:n
    %simulation of covariates
    x0 = 1;
    x1 = binornd(1,0.5);
    x2 = rand();
    X = [x0,x1,x2]';
    AllX(:,idx) = X;
end

%%%% simulated Betas
%simulate beta1 for the mean term
a = 2; b = 5;
beta0 = 20*s.^2;
%beta0 = a*zeros(1, m);
%beta0 = a*(-s+1);
beta2 = 20*sin(6*s) + a*s.^3 + a*ones(1, m);
beta1 = 15*cos(4*s).^4;
Beta_mu = [beta0;beta1;beta2];

figure(1), clf;
plot(s,beta0,'b','linewidth',3);
hold on;
plot(s,beta1,'g','linewidth',3);
hold on;
plot(s,beta2,'r','linewidth',3);
set(gca,'fontsize',22);
axis([0 1 -20 35])


%simulate beta2 for the variance term
%beta0 = ones(1,m);
beta0 = 0.2*(10 + 6*s.^3+ 3*cos(3*s));
a = 0.6;
beta2 =  a*s.^3 + a*ones(1, m);
beta3 = 0.8*s.^2;
Beta_sigma = [beta0;beta2;beta3];

figure(2), clf;
plot(s,beta0,'b','linewidth',3);
hold on;
plot(s,beta2,'g','linewidth',3);
hold on;
plot(s,beta3,'r','linewidth',3);
set(gca,'fontsize',22);
axis([0 1 0 5])

%%% parameters for simulation of noise
% senario 1, normal distribution, small noise
sigmae = 0.1;

%  senario 2, normal distribution, large noise
%sigmae = 0.4;

Mue = zeros(1,m);
Sigma_e = sigmae*eye(m);


Mu = zeros(1,m);
h = 0.8;
a = 0.6;
Sigma_gp = zeros(m,m);
for i=1:m
    for j=1:m
        %if i==j
        Sigma_gp(i,j) = a*exp(-((s(i)-s(j))/h)^2 );
        %end
    end
end


Sigma_e = Sigma_e/1;
Sigma_gp = Sigma_gp/1;
    

% parameters for coupola model
nu = 5/2;
Sigma = zeros(m,m,n);
h_min = s(2) - s(1);
alpha = [0.8,0.8,0.8];

%calculate the variance at each sample point
for i=1:n
    X = AllX(:,i);
    v_sigma = Beta_sigma'*X;
    for j=1:m
        Sigma(j,j,i) = (v_sigma(j))^2;
    end
    
    if(dispidx ==1)
        figure(10);clf;
        plot((v_sigma).^2);
        title('marginal variance');
        pause(0.01);
    end
end

%calculate the covariance
for idx = 1:n
    alpha_x = exp(alpha*AllX(:,idx));
    
    for i = 1:m
        for j=i+1:m
            besselj(nu,alpha_x*(j-i)*h_min);
            Sigma(i,j,idx) = sqrt(Sigma(i,i,idx))*sqrt(Sigma(j,j,idx))*(2^(1-nu)/gamma(nu))*(alpha_x*(j-i)*h_min)^nu*besselk(nu,alpha_x*(j-i)*h_min);
        end
    end
    
    Sigma(:,:,idx) = squeeze(Sigma(:,:,idx)) + squeeze(Sigma(:,:,idx))';
    
    
    for i=1:m
        Sigma(i,i,idx) = Sigma(i,i,idx)/2;
    end
    
    if(dispidx ==1)
        figure(10);clf;
        imagesc(Sigma(:,:,idx));
        colorbar;
        pause(0.1);
    end
    
    %sample from copula
    Sampled_Phi_Inv(idx,:) = mvnrnd(zeros(1,m),squeeze(Sigma(:,:,idx))); 
    for j=1:m
        Sampled_U(idx,j) = normcdf(Sampled_Phi_Inv(idx,j),0,sqrt(Sigma(j,j,idx)));
    end
    
end

%randomly sample U;
if(rand_indicator==1)
    Sampled_U = rand(size(Sampled_U));
end
%simulate marginal distribution


%plot the simulated Beta
taua = 0.01:0.01:0.99;
for idx_tau = 1:length(taua)
    tau = taua(idx_tau);
    Beta(:,:,idx_tau) = Beta_mu + Beta_sigma*norminv(tau,0,1);
    
    %plot the Beta
    if(dispidx ==1)
        figure(1),clf;
        plot(Beta(:,:,idx_tau)');
        pause(0.01);
    end
    
    %the marginal distribution of Y
    for i=1:n
        X = AllX(:,i);
        tmp_beta = squeeze(Beta(:,:,idx_tau));
        %add noise to the tau-th quantile
        
        QantileY(i,:,idx_tau) = tmp_beta'*X;
    end
end


% %calculate the variance at each sample point
% for i=1:n
%    v_sigma = Beta_sigma'*X; 
% end


% %the marginal distribution
% idx = 1;
% taua = 0.01:0.01:0.99;
% for idx_tau = 1:length(taua)
%     tau = taua(idx_tau);
%     for i=1:n
%         X = AllX(:,i);
%         muX = Beta_mu'*X;
%         stdX = sqrt(diag(Sigma(:,:,i))); %Beta_sigma'*X;
% 
%         QantileY(i,:,idx_tau) = muX + stdX.*norminv(tau,0,1);
%         
%     end
% end

%display the joint data
figure(1);clf;
displaytau = [50];
for i=1:length(displaytau)
   tau = displaytau(i);
   if(dispidx ==1)
       figure(1);hold on;
       plo t(s,squeeze( QantileY(:,:,tau)'),'linewidth',1.5);
       pause(0.5);
   end
end
axis([0 1 -10 40]);
set(gca,'fontsize',22);


%varify the std
for i=1:n
   simulatedY =  squeeze(QantileY(i,:,:));
   std_y = std(simulatedY');
   if(dispidx ==1)
       figure(1);clf;
       plot(std_y);
       hold on, plot(sqrt(diag(Sigma(:,:,i))),'-r','linewidth',2);
       pause(0.1);
   end
end

%the simulated data;
for i=1:n
    X = AllX(:,i);
    muX = Beta_mu'*X;
    stdX = Beta_sigma'*X;
    
    SimulatedY(i,:) = muX + stdX.*norminv(squeeze(Sampled_U(i,:))',0,1);
end


%display the simulated data
figure(2);clf;
plot(s,squeeze(SimulatedY'),'linewidth',1.5);
axis([0 1 -10 40]);
set(gca,'fontsize',22);

Ture_Sigma = Sigma;

%save simulated_1D_dataset1 SimulatedY QantileY Beta_mu Beta_sigma Sampled_U AllX m s n Ture_Sigma;


