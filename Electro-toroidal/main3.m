clc,clear all;
% % Define parameters:
% % Major and minor radii:
R_s = 0.4;
R_b = 1;
% % Parameters for compute non-dimensional presure P and electric load epsilon:
% P_tilde = 10;
% C_1 = ;
% C_2 = ;
H = 0.01;
% phi_0 = 1;
% beta = 1;
% % Global avariables:
global epsilon alpha gamma t_span rho0
% Pstart = (P_tilde * R_b) / (C_1 * H);
% epsilon = phi_0^2 / (C_1 * beta * H^2);
% % epsilon = 0;
% alpha = C_2/C_1;
% % aspect ratio
% gamma = R_s/R_b;
rho0 = 1.5;
rho_step = 0.01;
rho0_vec = [];
epsilon = 0.0;
alpha = 0.3;
gamma = R_s/R_b;
detai0 = 0.5;
P0 = 3;
%%%
% figure
% Initial guess:
% p_step = 0.1;
% epsilon_step = 0.01;
 % range of theta:
 n_elements = 500;
t_span = linspace(0,pi,n_elements);



