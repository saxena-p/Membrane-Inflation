function Sf = optfn(x)
% x = [detai0; Pi]
% For the ode45 calculation, y0 = [rho(0) rho'(0) eta(0) eta'(0)] = [rhoi0 0 0 detai0 ]
%     tic
    global drhoat0 h chixMU mus alphas r0 rrange rvals solution
%     global aig

    [rhoat1, etaat1, rvals, solution] = solnIVP(h, x(2), chixMU, mus, alphas, [r0*drhoat0; 0; x(1); 1e-8], rrange);
%     lambda10 = lambda20 = drhoat0;

    Sf = sqrt((rhoat1-1)^2 + etaat1^2);
%     clear diff
%     to
    
end