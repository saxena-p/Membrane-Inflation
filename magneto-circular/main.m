%% December 2017
%% Narravula Harshavardhan Reddy, GALCIT
%% Refer to the circular membrane paper for solution procedure and notations

clear
clc

global drhoat0 h chixMU mus alphas r0 rrange rvals solution

%% Material properties
mus = [1.4910 0.0029 -0.0236]; alphas = [1.3 5.0 -2.0]; 
%%we used \mu = % 0.414 MPa to arrive at these non-dimensional muk values; from Ogden 1972

% alpha = 0.00;
% mus = [1/(1-alpha) alpha/(1-alpha) 0]; alphas = [2 0 0];

chi = 2.5;
%% Dipole properties
MU = 0;
h = 1;  %non-dimensional height of the dipole

%% drhoat0 = lambda1 = lambda2 at r = 0;    % initial guess
drhoat0 = 2.0;

%% Initial guesses
Pig = 1.9038;    % internal pressure
eta0 = 1.0978;     % height of the deformed membrane at its pole

% Pig = 7;
% y0 = [detai0; rhoi0; 0];

chixMU = chi*MU;

r0 = 1e-6;
rrange = linspace(r0,1,1000);
% 
% *************************************
Sf = 1;
emax = 1e-10;

% tic
x0 = [eta0, Pig];

% while Sf > emax
% 
%     Sf1 = Sf;
%     tic
%     rhoi0 = rhoi0 + rhostep;

% % % % % % % % % % % % % % % % % % %     
%     NOTE: tolerances here should match with those of ode solver in
%     solnIVP.m or it doesn't make sense!
     options = optimset('FunValCheck', 'on', 'MaxFunEvals', 3500, 'TolFun', 1e-13, 'TolX', 1e-05);
%     options = optimset('FunValCheck', 'on', 'MaxFunEvals', 10000, 'TolFun', 1e-13, 'TolX', 1e-06);
%                      'PlotFcn', @optimplotfval,'PlotFcn', @optimplotx, 
%     [x, fval, exitflag, output] = fminsearch(@optfn, x0, options);
    [x, fval, exitflag, output] = fminsearchbnd(@optfn, x0, [0, 0], [], options);
%     options = optimset('MaxFunctionEvaluations', 1000);
%     [x, fval, exitflag, output] = fmincon(@optfn, x0, [], [], [], [], [0, 0], [Inf, Inf],[],options);
    Sf = fval;
%     toc
    
%     if Sf > Sf1
%        disp('Sf increases here!');
%        flag = 1;
%        break
%     else
%         x0 = x;     %% if Sf increases, this will avoid taking x corresponding to the higher Sf value
%         flag = 0;
%     end
% end

% toc

[drhoat1, etaat1, rvals, soln] = solnIVP(h, x(2), chixMU, mus, alphas, [r0*drhoat0; 0; x(1); 0], rrange);     %% if Sf increases, this will avoid taking x corresponding to the higher Sf value
Sf_check = sqrt(drhoat1^2 + etaat1^2);
lambda2 = soln(:,1)./rvals;
lambda1 = sqrt((soln(:,2)+lambda2).^2 + soln(:,4).^2);
% a = a;

lambda22 = solution(:,1)./rvals;
lambda11 = sqrt((solution(:,2)+lambda22).^2 + solution(:,4).^2);


% figure('PaperPositionMode','auto');
% set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'FontWeight', 'normal', 'FontName', 'Times');

out = [drhoat0 x(1) x(2)]
Sf

hold on
plot(solution(:,1),solution(:,3),'LineWidth',0.75)
xlabel('$\varrho$','Interpreter','LaTeX','FontSize', 16, 'FontWeight', 'normal', 'FontName', 'Times');
ylabel('$\eta$','Interpreter','LaTeX','FontSize',16, 'FontWeight', 'normal', 'FontName', 'Times');

%% Second variation:

[u, rs, detU, P11, P22, detPs, posdefs] = second_variation(rvals, solution, lambda22, h, x(2), chixMU, mus, alphas);

nansu = isnan(u);
if any(nansu(:)==1) == 1
   disp('Warning: NANs in u!') 
end

suff_check1 = 0;
for i = 1:length(u)
   if u(i,1)==0 && u(i,2)==0 && u(i,3)==0 && u(i,4)==0
      suff_check1 = 1;
   end
end
suff_check1
suff_check2 = any(detU(2:end)==0)
ness_check = any(posdefs>0)

%     rhoval = solution(:,1); drhoval = lambda22 + solution(:,2);
%     etaval = solution(:,3); detaval = solution(:,4);
%     
%     rho_interp = griddedInterpolant(rvals, rhoval, 'pchip'); eta_interp = griddedInterpolant(rvals, etaval, 'pchip');
%     drho_interp = griddedInterpolant(rvals, drhoval, 'pchip'); deta_interp = griddedInterpolant(rvals, detaval, 'pchip');
%     
%     rs = (r0:1e-5:1);
%     etas = zeros(size(rs));
%     
%     for i = 1:length(rs)
%        etas(i) = eta_interp(rs(i)); 
%     end
%     
%     figure;
%     plot(rvals, etaval, '-b', rs, etas, '-rx')