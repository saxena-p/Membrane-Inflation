% Rewriting a clearer and better commented version of code originally
% written by Harsha.

clear
% clc

% define global variables that are used throughout
global drhoat0 h chixMU mus alphas r0 rrange rvals solution

% Material properties (Ogden elastic parameters)
mus = [1.4910 0.0029 -0.0236]; 
alphas = [1.3 5.0 -2.0]; 

% Coupling parameter chi
chi = 2.5;

% Dipole properties
MU = 0;
h = 1;  %non-dimensional height of the dipole

%% Fixing lambda_2 and lambda_1
drhoat0 = 2.0;  % lambda value at the centre point

%% initial guesses
Pig = 1.8038;   % internal pressure
eta0 = 1.0978;  % height of the deformed membrane at its pole


chixMU = chi*MU;

% instead of starting at 0, we start integration at r0.
r0 = 1e-6;
rrange = linspace(r0,1,1000);

% Do the optimization
x0 = [eta0, Pig];
options = optimset('FunValCheck', 'on', 'MaxFunEvals', 3500, 'TolFun', 1e-13, 'TolX', 1e-05);
[x, fval, exitflag, output] = fminsearchbnd(@optfn, x0, [0, 0], [], options);
Sf = fval;  % Variable to store the final optimized function value
disp(['The optimized function value is ', num2str(Sf,5)]);
disp(['Optimized eta_0 is ', num2str(x(1))]);
disp(['Optimized pressure is ', num2str(x(2))]);

% Get solution for the optimized profile
[drhoat1, etaat1, rvals, soln] = solnIVP(h, x(2), chixMU, mus, alphas, [r0*drhoat0; 0; x(1); 0], rrange);     %% if Sf increases, this will avoid taking x corresponding to the higher Sf value
Sf_check = sqrt(drhoat1^2 + etaat1^2);
lambda2 = soln(:,1)./rvals;
lambda1 = sqrt((soln(:,2)+lambda2).^2 + soln(:,4).^2);

% Plot the membrane profile
figure;
plot(solution(:,1),solution(:,3),'LineWidth',0.75)
xlabel('$\varrho$','Interpreter','LaTeX','FontSize', 16, 'FontWeight', 'normal', 'FontName', 'Times');
ylabel('$\eta$','Interpreter','LaTeX','FontSize',16, 'FontWeight', 'normal', 'FontName', 'Times');
grid on

% Plot the lambda values.
figure;
plot(rvals, lambda1,'LineWidth', 2)
hold on
plot(rvals, lambda2,'LineWidth', 2)
legend('\lambda_1', '\lambda_2')
xlabel('radius')
ylabel('Stretch values')
grid on

% Write the data on files
