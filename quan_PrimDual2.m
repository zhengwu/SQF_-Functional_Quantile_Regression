function [bta, dd_m, Yhat, gacv, bceof,lagrange_value] = quan_PrimDual2(Y, X, Sig, tau, lam)
% Quantile function-on-scalar regression
% Primal dual algorithm
% Yhat --- fitted value
% gacv --- generalized approximate CV
% rho_tau(Y-hat Y)/(mn - df)

% lagrange_value - output

[n, m] = size(Y);
[p, n1] = size(X);
if n1~=n
    disp('dimension error!')
end
nm = n*m;

%lb = -(1-tau)*ones(nm, 1); ub = tau*ones(nm, 1);
lb = -(1-tau); ub = tau;
XX = X'*X;       % size n X n
%Q = kron(Sig, XX);    % size nm X nm   %% waste a lot of time!!!
%Xtilde = kron(ones(m, 1), X')';  % p X nm
Yvec = Y(:);

mu0 = (X*X')\ (X* mean(Y, 2));
%mu0 = zeros(p, 1);
J = 3;
k=0;

while k <= J
    
    %k
    k = k+1;
    
    Xt_mu0 = X'*mu0;
    Xtilde_mu0 = repmat(Xt_mu0, m, 1);
    
    %bb = lam*(Yvec-Xtilde'*mu0);
    bb = lam*( Yvec - Xtilde_mu0  );
    
    %dd = quad_box(Q, bb, lb, ub);
    dd = quad_box2(Sig, XX, bb, lb, ub);
    dd_id = find( (dd > lb).*(dd<ub)  ); 
    
    n_dd = length(dd_id);
    if isempty(dd_id) || n_dd < p 
        %disp('please change the initial value!');
        %break
        dd_id = unidrnd(nm, p+1, 1);
        n_dd = p+1;
    end
    
    dd_fset = floor((dd_id-0.5)/n);
    dd_set1 = dd_id - dd_fset*n;  %%% position of row
    dd_set2 = dd_fset + 1;  %%% position of col

    Xdd = X(:, dd_set1);  Xdd2 = Xdd*Xdd';
    Ydd = Y(dd_id);
    
    dd_m = reshape(dd, n, m);
    %%%Sig_XX_dd = XX(dd_set1, :)*dd_m*Sig(dd_set2, : )';
    Sig_XX_dd = zeros(n_dd, 1);
    for j = 1:n_dd
        %%%%%Sig_XX_dd(j, :)  = kron( Sig(dd_set2(j), : ), XX(dd_set1(j), :) );
        Sig_XX_dd(j) = XX(dd_set1(j), :)*dd_m*Sig(dd_set2(j), : )';
    end
    
    %mu0 = (Xdd2) \ ( Xdd *(  Ydd - (Sig_XX_dd*dd)/(2*lam) ) );
    %mu0 = (Xdd2) \ ( Xdd *(  Ydd - (Sig_XX_dd)/(2*lam) ) );
    
    mu0 = (Xdd2) \ ( Xdd *(  Ydd - (Sig_XX_dd)/(lam) ) ); % fixed by Zhengwu
end

bta = zeros(p, m);
bcoef = zeros(p,m);

for k = 1:p
    %for j = 1:m
    %    bta(k, j) = mu0(k)  + (1/lam)* dd'*kron( Sig(:, j), X(k, :)'  );
    %end
    bta(k, :) = ( mu0(k)*ones(m, 1) + (1/lam)* (kron(Sig, X(k, :)'  )'*dd ) )';
    bceof(k,:) =  (1/lam)* (kron(eye(size(Sig)), X(k, :)'  )'*dd );
end
%plot(bta'); hold on; plot(Beta', 'r--'); hold off;

Yhat = X'*bta;
Yhat_vec = Yhat(:);
Ydiff = Yvec - Yhat_vec;
Ydiff_id = (Ydiff<0);
gacv =  sum( Ydiff.*(tau - Ydiff_id) ) /(nm - n_dd);

lagrange_value.mu0 = mu0;
lagrange_value.dd = dd;


