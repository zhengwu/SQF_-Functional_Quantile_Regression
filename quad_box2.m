function x1 = quad_box2(Sig, XX, b, l, u)
%%%
% min 1/2 x'Gx - b'x
% G is a kron(Sig, XX)
% l \le x \le u

%L = sum( diag(G) );
L = trace(Sig)*trace(XX);

%[n1, n2] = size(G);
%if n1~=n2
%    disp('dimension error!');
%end
[m, ~] = size(Sig);
[n, ~] = size(XX);


x0 = l;
x1 = b/L;

z1 = x1; 
k =1; 
t1 =1;

J = 200;

while max(abs(x1-x0))>1e-4 && k<=J  

    
    k = k+1;
    
    x0 = x1;
    
    %tic
    z1_m = reshape(z1, n, m);
    Gz1_m = XX*z1_m*Sig';
    Gz1 = Gz1_m(:);
    x1_nobox = 1/L*(L*z1 + b- Gz1);    
    x1 = min(u, max(l, x1_nobox));
    %toc
    
    t2 = ( 1+sqrt(1 + 4*t1^2))/2;
    z2 = x1 + (t1-1)/t2*( x1 - x0);
    
 
    z1 = z2;
    t1 = t2;
    
    %plot(k, 1/2*x1'*G*x1 - x1'*b, 'x'); hold on;
    %plot(k, max(abs(x1-x0)), 'x'); hold on;
    
end




