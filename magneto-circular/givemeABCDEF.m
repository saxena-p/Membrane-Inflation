function [A, B, C, D, E, F] = givemeABCDEF(y, lambda2, r, h, P, MU, mus, alphas)

%% Refer to the section Solution Procedure in the paper
% % dWdl1 is the result after dividing by \mu, shear modulus

% % % rho = y(1); drho = y(2); eta = y(3); deta = y(4);
% % % u = rho; v = drho; w = eta; x = deta;

% % % l1 = sqrt(drho^2 + deta^2); l2 = rho/r; l3 = 1/(l1*l2);

u = y(1); v = y(2); eta = y(3); w = y(4);

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

a = zeros(1,3); b = zeros(1,3); c = zeros(1,3);
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

% % % L11 = v/l1; L12 = x/l1; L22 = v/r - u/r^2; L31 = -v/(l2*l1^3); L32 = -x/(l2*l1^3); L33 = -L22/(l1*l2^2);
% % % 
% % % W11 = 0; W12 = 0; W13 = 0;
% % % for k = 1:3
% % %    
% % %     W11 = W11 + mus(k)* ( (alphas(k)-1)*(l1^(alphas(k)-2))*L11 - (alphas(k)/l1)*(l3^(alphas(k)-1))*L31 + (l3^alphas(k)/ l1^2)*L11 );
% % %     W12 = W12 + mus(k)* ( (alphas(k)-1)*(l1^(alphas(k)-2))*L12 - (alphas(k)/l1)*(l3^(alphas(k)-1))*L32 + (l3^alphas(k)/ l1^2)*L12 );
% % %     W13 = W13 + mus(k)* ( -(alphas(k)/l1)*(l3^(alphas(k)-1))*L33 );
% % %     
% % % end
% % % 
% % % A = (r/l1)*dWdl1 - (r*v/l1^2)*dWdl1*L11 + (r*v/l1)*W11;
% % % B = -(r*v/l1^2)*dWdl1*L12 + (r*v/l1)*W12;
% % % C1 = (v/l1)*dWdl1 + (r*v/l1)*W13;
% % % length1 = rho^2 + 5*(eta-h)^2; length2 = rho^2 + (eta-h)^2;
% % % C = C1 + dWdl2 + 3*MU*r*rho*length1/length2^5 + P*u*x;
% % % D = -(r*x/l1^2)*dWdl1*L11 + (r*x/l1)*W11;
% % % E = (r/l1)*dWdl1 - (r*x/l1^2)*dWdl1*L12 + (r*x/l1)*W12;
% % % F1 = (x/l1)*dWdl1 + (r*x/l1)*W13;
% % % F = F1 + 12*MU*r*(eta-h)^3/length2^5 - P*u*v;

end