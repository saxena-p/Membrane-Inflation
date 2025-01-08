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
rho_step = 0.01;
rho0_vec = [];
epsilon = 0.04;
alpha = 0.3;
gamma = R_s/R_b;
file_check = fopen('initial_guess.txt','r');
pm = str2num(fgetl(file_check));
rho0 = 1.5;
fclose(file_check);
if(pm(1) == gamma && pm(2)==epsilon && pm(3) == alpha)
    file_read = fopen('initial_guess.txt','r');
    fileID = fopen('initial_guess.txt','a+');
    line = fgetl(file_read);
    line = fgetl(file_read);
    data = [];
    line = fgetl(file_read);
    while(line ~= -1)
        data = [data; str2num(line)];
        line = fgetl(file_read);
    end
    for i = 1:size(data,1)
        if rho0 == data(i,1)
            P0 = data(i,2);
            detai0 = data(i,3);
            break
        end
    end
else
    fileID = fopen('initial_guess.txt','w');
    fprintf(fileID,'%3.3f %3.3f %3.3f\n',[gamma;epsilon;alpha]);
    fprintf(fileID,'rho0 P0 deta0\n');
    rho0 = 1.5;
    detai0 = 0.5;
    P0 = 3;
    % rho0 = 4.828;
    % detai0 = 3.058;
    % P0 = 6.0720;
end
%%%
% figure
% Initial guess:
% p_step = 0.1;
% epsilon_step = 0.01;
 % range of theta:
 n_elements = 500;
t_span = linspace(0,pi,n_elements+1);
lambda1vec0 = [];
lambda2vec0 = [];
lambda1vechpi = [];
lambda2vechpi = [];
lambda1vecpi = [];
lambda2vecpi = [];
Pvec=[];
Vvec=[];
stress2vec = [];
warning off;
for i = 1:500
%  P= P+p_step
%  epsilon = epsilon + epsilon_step
rho0 = rho0+rho_step
for j = 1:size(data,1)
        if rho0 == data(j,1)
            P0 = data(j,2);
            detai0 = data(j,3);
            break
        end
end
%  Pvec = [Pvec P];
 y0 = [rho0;0.0;0;detai0];
 
 x0=[P0,detai0];

 options = optimset('FunValCheck', 'on', 'MaxFunEvals', 10000, 'TolFun', 1e-12, 'TolX', 1e-06);

 [x, fval, exitflag, output] = fminsearchbnd(@opt_fun_2, x0, [0,0], [], options);
 Sf = fval;
 
 [drhoatpi, etaatpi, t, y] = IVP_solver(gamma, alpha, x(1), epsilon, [rho0;0;0;x(2)], t_span);
 Pvec = [Pvec x(1)];
 figure(7)
 plot(y(:,1),y(:,3), 'linewidth', 0.5);
 xlim([0,6]);
 ylim([0,6]);
 hold on;
 lambda1 = sqrt(y(:,2).^2 + y(:,4).^2)/gamma;
 lambda2 = y(:,1)./(1+gamma*cos(t));
 lambda3 = 1./(lambda1.*lambda2);
 s11 = - (x(1) * H)/R_b + 2 .* lambda1 .* lambda1 +...
     2*alpha .* lambda1 .* lambda1 .* (lambda2.*lambda2 + 1./(lambda1 .* lambda1 .* lambda2 .* lambda2))...
     -2./(lambda1.*lambda1.*lambda2.*lambda2) - 2*alpha.*(1./(lambda1.*lambda1) + 1./(lambda2.*lambda2))...
     -0.5 * epsilon .* lambda1 .* lambda1 .* lambda2 .* lambda2;
 s22 = - (x(1) * H)/R_b + 2 .* lambda2 .* lambda2 +...
     2*alpha .* lambda2 .* lambda2 .* (lambda1.*lambda1 + 1./(lambda1 .* lambda1 .* lambda2 .* lambda2))...
     -2./(lambda1.*lambda1.*lambda2.*lambda2) - 2*alpha.*(1./(lambda1.*lambda1) + 1./(lambda2.*lambda2))...
     -0.5 * epsilon .* lambda1 .* lambda1 .* lambda2 .* lambda2
 for i = 1:size(s22,1)
     if s22(i) < 0
         theta_limit = i*pi/500
         break
     end
 end
 lambda1vec0 = [lambda1vec0 lambda1(1)];
 lambda2vec0 = [lambda2vec0 lambda2(1)];
 lambda1vechpi = [lambda1vechpi lambda1(n_elements/2)];
 lambda2vechpi = [lambda2vechpi lambda2(n_elements/2)];
 lambda1vecpi = [lambda1vecpi lambda1(end)];
 lambda2vecpi = [lambda2vecpi lambda2(end)];
 P0 = x(1);
 
 exisit = 0;
 for i = 1:size(data,1)
     if abs(rho0 - data(i,1)) < 10e-6
         exisit = 1;
         break
     end
 end
 if exisit == 0
    fprintf(fileID,'%3.5f %3.5f %3.5f\n',[rho0; P0; detai0]);
 end

 stress1 = s11(end-1);
 stress2 = s22(end-1);
 stress2vec = [stress2vec stress2];
 rho0_vec = [rho0_vec rho0];
 x
 volume = 0;
 for ip = 2:n_elements+1
     rho_a = y(ip-1,1);
     eta_a = y(ip-1,3);
     rho_b = y(ip,1);
     eta_b = y(ip,3);
     volume = volume + pi*rho_a*(eta_a+eta_b)*(rho_a-rho_b);
 end
 volume
 Vvec = [Vvec volume];
end
fclose(fileID);
figure(1)
plot(lambda1vec0,Pvec);
hold on;
figure(2)
plot(lambda1vechpi,Pvec);
hold on;
figure(3)
plot(lambda1vecpi,Pvec);
hold on;
figure(4)
plot(lambda2vec0,Pvec);
hold on;
figure(5)
plot(lambda2vechpi,Pvec);
hold on;
figure(6)
plot(lambda2vecpi,Pvec);
hold on;
figure(8)
plot(rho0_vec,stress2vec);
hold on;