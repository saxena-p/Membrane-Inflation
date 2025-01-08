function [u, rs, detU, P11, P22, detPs, posdefs] = second_variation(rval, solution, lambda2, h, P, chixMU, mus, alphas)    


    rhoval = solution(:,1); drhoval = lambda2 + solution(:,2);
    etaval = solution(:,3); detaval = solution(:,4);
    
%     [AA, BB, CC, DD, EE, FF] = ode45essentials(solution, lambda2, rval, h, P, chixMU, mus, alphas);
%     
%     d2rhoval = zeros(size(rhoval)); d2etaval = zeros(size(etaval));
%     for i = 1:length(d2rhoval)
%        
%         Amat = [AA(i) BB(i); DD(i) EE(i)]; Emat = [CC(i); FF(i)];
%         dvdw = Amat\Emat;
%         d2rhoval(i) = dvdw(1); d2etaval(i) = dvdw(2);
%         
%     end
    
    d2rhoval = zeros(size(rhoval)); d2etaval = zeros(size(etaval));
    dr = rval(2) - rval(1);
    nodes = length(d2rhoval);
    d2rhoval(1) = (drhoval(2) - drhoval(1))/(dr); d2rhoval(nodes) = (drhoval(nodes) - drhoval(nodes-1))/dr;
    d2etaval(1) = (detaval(2) - detaval(1))/(dr); d2etaval(nodes) = (detaval(nodes) - detaval(nodes-1))/dr;
    for i = 2: nodes-1

        d2rhoval(i) = (drhoval(i+1) - drhoval(i-1))/(2*dr);
        d2etaval(i) = (detaval(i+1) - detaval(i-1))/(2*dr);

    end
    
    d3rhoval = zeros(size(rhoval)); d3etaval = zeros(size(etaval));
    
    dr = rval(2) - rval(1);
    nodes = length(d3rhoval);
    d3rhoval(1) = (d2rhoval(2) - d2rhoval(1))/(dr); d3rhoval(nodes) = (d2rhoval(nodes) - d2rhoval(nodes-1))/dr;
    d3etaval(1) = (d2etaval(2) - d2etaval(1))/(dr); d3etaval(nodes) = (d2etaval(nodes) - d2etaval(nodes-1))/dr;
    for i = 2: nodes-1

        d3rhoval(i) = (d2rhoval(i+1) - d2rhoval(i-1))/(2*dr);
        d3etaval(i) = (d2etaval(i+1) - d2etaval(i-1))/(2*dr);

    end
    
    rho_interp = griddedInterpolant(rval, rhoval, 'pchip'); eta_interp = griddedInterpolant(rval, etaval, 'pchip');
    drho_interp = griddedInterpolant(rval, drhoval, 'pchip'); deta_interp = griddedInterpolant(rval, detaval, 'pchip');
    d2rho_interp = griddedInterpolant(rval, d2rhoval, 'pchip'); d2eta_interp = griddedInterpolant(rval, d2etaval, 'pchip');
    d3rho_interp = griddedInterpolant(rval, d3rhoval, 'pchip'); d3eta_interp = griddedInterpolant(rval, d3etaval, 'pchip');
    
    options = odeset('stats','off','Mass',@mass,'RelTol',1e-5,'AbsTol',1e-10);
    
    u0 = [0 1 0 0 0 0 0 1];
%     trange = linspace(0, pi, 200);
    [rs, u] = ode45(@secondvar_IVP,rval,u0,options);
    
    detU = u(:,1).*u(:,7) - u(:,3).*u(:,5);
    
    [P1, P2, P3, P4, ~, ~, ~, ~, ~, ~, ~, ~] = PdPQ([rhoval, drhoval, d2rhoval, d3rhoval, etaval, detaval, d2etaval, d3etaval], rval, h, P, chixMU, mus, alphas);
    P11 = P1; P22 = P4;
% %     Matrix is positive definite if the second out argument of chol is zero,
% not so if it is a positive integer
    
    detPs = P1.*P4 - P2.*P3;
    posdefs = ones(length(rhoval), 1);
    
    for i = length(posdefs)
        [~, posdefs(i)] = chol([P1(i,1) P2(i,1); P3(i,1) P4(i,1)]);
    end
    
    function dydt = secondvar_IVP(r,u)

% %         r
        rho = rho_interp(r);
        eta = eta_interp(r); drho = drho_interp(r); deta = deta_interp(r);
        d2rho = d2rho_interp(r); d2eta = d2eta_interp(r);
        d3rho = d3rho_interp(r); d3eta = d3eta_interp(r);

        [~, ~, ~, ~, dP1, dP2, dP3, dP4, Q1, Q2, Q3, Q4] = PdPQ([rho, drho, d2rho, d3rho, eta, deta, d2eta, d3eta], r, h, P, chixMU, mus, alphas);
% u
        dydt = [u(2); dP1*u(2) + dP2*u(6) - Q1*u(1) - Q2*u(5); u(4); dP1*u(4) + dP2*u(8) - Q1*u(3) - Q2*u(7); ...
            u(6); dP3*u(2) + dP4*u(6) - Q3*u(1) - Q4*u(5); u(8); dP3*u(4) + dP4*u(8) - Q3*u(3) - Q4*u(7)];

    end

    function M = mass(r,~)

        rho = rho_interp(r);
        eta = eta_interp(r); drho = drho_interp(r); deta = deta_interp(r);
        d2rho = d2rho_interp(r); d2eta = d2eta_interp(r);
        d3rho = d3rho_interp(r); d3eta = d3eta_interp(r);

        [P1, P2, P3, P4, ~, ~, ~, ~, ~, ~, ~, ~] = PdPQ([rho, drho, d2rho, d3rho, eta, deta, d2eta, d3eta], r, h, P, chixMU, mus, alphas);
        
        M = eye(8,8);
        M(2,2) = -P1; M(2,6) = -P2; M(4,4) = -P1; M(4,8) = -P2;
        M(6,2) = -P3; M(6,6) = -P4; M(8,4) = -P3; M(8,8) = -P4;

    end
    
end
    
    