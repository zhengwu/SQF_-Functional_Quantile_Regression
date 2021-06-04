function f = Qtheta(theta)

global Param;

X = Param.X;
nu = Param.nu;
two_gamma_hat_h = Param.two_gamma_hat_h;
V = Param.V;
h = Param.h;



[Nh,n] = size(two_gamma_hat_h);


p = ceil(nu)-1;
if(nu-p ~= 0.5)
    warning('nu must be i/2...');
    return;
end

gx = zeros(Nh,n);

if(length(h) == 1)
    h_min = h;
else
    h_min = h(2) - h(1);
end;

for ih = 1:Nh
    for i=1:n
        curr_x = X(:,i);        
        gx(ih,i) = two_gamma_hat_h(ih,i) - 2 + 2*(2^(1-nu)/gamma(nu))*(exp(curr_x'*theta)*ih*h_min)^nu*besselk(nu,exp(curr_x'*theta)*ih*h_min);
    end
end

f = 0;

for i = 1:n
   f = f + gx(:,i)'*V*gx(:,i); 
end


