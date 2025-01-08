function Sf = opt_fun(x)
    global P epsilon alpha gamma t_span
    [drhoatpi, etaatpi, t, y] = IVP_solver(gamma, alpha, P, epsilon, [x(1),0,0,x(2)], t_span);
    Sf = sqrt(drhoatpi^2 + etaatpi^2);
%     toc
end

