
clear
clc

syms u v w eta lambda2 h P MU r mu1 mu2 mu3 alpha1 alpha2 alpha3

mus = [mu1 mu2 mu3]; alphas = [alpha1 alpha2 alpha3];

up = lambda2 + v;
lambda1 = sqrt(up^2 + w^2);
rho = u;

lambda3 = 1/(lambda1*lambda2);

dWdl1 = 0; dWdl2 = 0;
for k = 1:3
   
    dWdl1 = dWdl1 + mus(k)*( lambda1^(alphas(k) - 1) - lambda3^alphas(k)/lambda1 );
    dWdl2 = dWdl2 + mus(k)*( lambda2^(alphas(k) - 1) - lambda3^alphas(k)/lambda2 );
    
end
clear k

a = sym('a', [1 3]); b = sym('b', [1 3]); c = sym('c', [1 3]);
for k = 1:3
    
    a(k) = mus(k)*(alphas(k) - 1)*lambda1^( alphas(k)-2 );
    b(k) = mus(k)*(1/lambda2)*alphas(k)*lambda3^( alphas(k)-1 );
    c(k) = mus(k)*(1/lambda2^2)*lambda3^alphas(k);
    
end
clear k

d = a/lambda1 + b/(lambda2*lambda1^3); e = r*(w/lambda1)*a + (r*w/(lambda2*lambda1^3))*b;

W11 = 0; W12 = 0; W13 = 0;
for k = 1:3
    
    W11 = W11 + r*up*d(k);
    W12 = W12 + e(k);
    W13 = W13 + up*v*d(k) + v*( b(k)/(lambda1*lambda2^2) + c(k) );
    
end
clear k

WW11 = dWdl1*(r/lambda1^2)*(lambda1 - up^2/lambda1) + (up/lambda1)*W11;
WW12 = dWdl1*(1/lambda1^2)*(-r*w*up/lambda1) + (up/lambda1)*W12;
WW13 = (up/lambda1)*dWdl1 + dWdl1*(1/lambda1^2)*v*(lambda1 - up^2/lambda1) + (up/lambda1)*W13;

WW21 = -r*w*dWdl1*(up/lambda1^3) + (w/lambda1)*W11;
WW22 = (r/lambda1)*dWdl1 - (r*w^2/lambda1^3)*dWdl1 + (w/lambda1)*W12;
WW23 = (w/lambda1)*dWdl1 - (w/lambda1^3)*dWdl1*up*v + (w/lambda1)*W13;

A = WW11; B = WW12;
length1 = rho^2 + 5*(eta-h)^2; length2 = rho^2 + (eta-h)^2;
C = -WW13 + dWdl2 + 3*MU*r*rho*length1/length2^5 + P*u*w;

D = WW21; E = WW22;
F = -WW23 + 12*MU*r*(eta-h)^3/length2^5 - P*u*up;

j = matlabFunction(A,B,C,D,E,F,...
    'file', 'ode45essentials.m', 'vars', {[u, v, eta, w] lambda2 r h P MU [mu1 mu2 mu3] [alpha1 alpha2 alpha3]},...
    'outputs', {'A','B','C','D','E','F'});