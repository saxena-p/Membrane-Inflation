clc,clear all;
% % Define parameters:

% % Parameters for compute non-dimensional presure P and electric load epsilon:
% P_tilde = 10;

% phi_0 = 1;
% beta = 1;
% % Global avariables:
global epsilon alpha gamma t_span rho0 H R_b
% thickness
H = 0.01;
% % Major and minor radii:
R_s = 0.4;
R_b = 1;

% electrical load : epsilon = phi_0^2 / (C_1 * beta * H^2);
epsilon = 0.04;
% alpha = C_2/C_1;
alpha = 0.3;
% aspect ratio of radii
gamma = R_s/R_b;

% another option
% rho0 = 4.49;
% detai0 = 2.3873;
% P0 = 5.8021;

% initial guess
rho0 = 4.828;
detai0 = 3.058;
P0 = 6.0720;
% step size of rho
rho_step = 0.005;


 % range of theta:
n_elements = 500;
t_span = linspace(0,pi,n_elements);

 % initialise some vectors
rho0_vec = [];
lambda1vec0 = [];
lambda2vec0 = [];
lambda1vechpi = [];
lambda2vechpi = [];
lambda1vecpi = [];
lambda2vecpi = [];
Pvec=[];
stress2vec = [];
warning off;

for i = 1:5
%  P= P+p_step
%  epsilon = epsilon + epsilon_step
rho0 = rho0+rho_step
%  Pvec = [Pvec P];
 y0 = [rho0;0.0;0;detai0];
 
 x0=[P0,detai0];

 options = optimset('FunValCheck', 'on', 'MaxFunEvals', 10000, 'TolFun', 1e-12, 'TolX', 1e-06);
 
 % we first solve the problem without energy relaxation to find a
 % reasonable initial guess
 
%  [x, fval, exitflag, output] = fminsearchbnd(@opt_fun_ER, x0, [0,0], [], options);
 [x, fval, exitflag, output] = fminsearchbnd(@opt_fun_2, x0, [0,0], [], options);

 Sf = fval;
 
 [drhoatpi, etaatpi, t, y] = IVP_solver(gamma, alpha, x(1), epsilon, [rho0;0;0;x(2)], t_span, H, R_b);

 Pvec = [Pvec x(1)];
 figure(1)
 plot(y(:,1),y(:,3), 'linewidth', 0.5);
 xlim([0,6]);
 ylim([0,6]);
 hold on;
 
 % here we use the results as the initial guess for solving the problem
 % with energy relaxation in the region will S_22 is negative.
   
 P0 = x(1);
 detai0 = y(1,4);
 [x2, fval, exitflag, output] = fminsearchbnd(@opt_fun_ER, x0, [0,0], [], options);
  Sf = fval;

 x0=[P0,detai0];
 [drhoatpi2, etaatpi2, t2, y2] = IVP_solver_ER(gamma, alpha, x2(1), epsilon, [rho0;0;0;x2(2)], t_span, H, R_b);
 
 figure(2)
 plot(y2(:,1),y2(:,3), 'linewidth', 0.5);
 xlim([0,6]);
 ylim([0,6]);
 hold on;
%  
 lambda1 = sqrt(y(:,2).^2 + y(:,4).^2)/gamma;
 lambda2 = y(:,1)./(1+gamma*cos(t));
 
  P0 = x2(1);
 detai0 = y2(1,4);

end
figure(3)
 plot(y(:,1),y(:,3),'-b', 'linewidth', 0.5);
 hold on;
 plot(y2(:,1),y2(:,3),'-r', 'linewidth', 0.5);
 hold on;
 xlim([0,6]);
 ylim([0,6]);