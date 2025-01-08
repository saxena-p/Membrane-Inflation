function Sf = opt_fun_ER(x)
    global epsilon alpha gamma t_span rho0 H R_b
    P = x(1);
    [drhoatpi, etaatpi, t, y] = IVP_solver_ER(gamma, alpha, P, epsilon, [rho0,0,0,x(2)], t_span, H, R_b);
    Sf = sqrt(drhoatpi^2 + etaatpi^2);
%     toc
end

