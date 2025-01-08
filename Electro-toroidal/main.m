clc,clear all;
% % Define parameters:
% % Major and minor radii:
% R_s = 0.2;
% R_b = 1;
% % Parameters for compute non-dimensional presure P and electric load epsilon:
% P_tilde = 10;
% C_1 = 100;
% C_2 = 10;
% H = 0.01;
% phi_0 = 1;
% beta = 1;
% % Global avariables:
global P epsilon alpha gamma t_span
% Pstart = (P_tilde * R_b) / (C_1 * H);
% epsilon = phi_0^2 / (C_1 * beta * H^2);
% % epsilon = 0;
% alpha = C_2/C_1;
% % aspect ratio
% gamma = R_s/R_b;

%%%
P = 4.9;
   
epsilon = 0.;
alpha = 0.2;
gamma = 0.2;
%%%
% Initial guess:
detai0 = 0.2;
rhoi0 = 1.4;
p_step = 0.1;
epsilon_step = 0.01;
 % range of theta:
t_span = linspace(0,pi,60);
Pvec = [];
lambda1vec0 = [];
lambda2vec0 = [];
lambda1vechpi = [];
lambda2vechpi = [];
lambda1vecpi = [];
lambda2vecpi = [];

figure

for i = 1:1
 P= P+p_step
%  epsilon = epsilon + epsilon_step
 Pvec = [Pvec P];
 y0 = [rhoi0;0.0;0;detai0];
 
 x0=[y0(1),y0(4)];
 
 options = optimset('FunValCheck', 'on', 'MaxFunEvals', 10000, 'TolFun', 1e-14, 'TolX', 1e-06);
 
 [x, fval, exitflag, output] = fminsearchbnd(@opt_fun, x0, [0,0], [ ], options);
 Sf = fval
 
 [drhoatpi, etaatpi, t, y] = IVP_solver(gamma, alpha, P, epsilon, [x(1);0;0;x(2)], t_span);
 lambda1 = sqrt(y(:,2).^2 + y(:,4).^2)/gamma;
 lambda2 = y(:,1)./(1+gamma*cos(t));
 lambda3 = 1./(lambda1.*lambda2);
 lambda1vec0 = [lambda1vec0 lambda1(1)];
 lambda2vec0 = [lambda2vec0 lambda2(1)];
 lambda1vechpi = [lambda1vechpi lambda1(30)];
 lambda2vechpi = [lambda2vechpi lambda2(30)];
 lambda1vecpi = [lambda1vecpi lambda1(end)];
 lambda2vecpi = [lambda2vecpi lambda2(end)];
 
 plot(y(:,1),y(:,3), 'linewidth', 0.5);
 hold on;
 rhoi0 = y(1,1);
 detai0 = y(1,4);
end