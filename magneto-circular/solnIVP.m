function [rhoat1, etaat1, rvals, sol] = solnIVP(h, P, chixMU, mus, alphas, y0, rrange)

% global drhoat0
    % refer to the notebook for the derivation. it's simple
    
%     options = odeset('Mass',@mass,'MStateDependence','strong');
%     NOTE: tolerances here should match with those of fminsearch in
%     main.m or it doesn't make sense!

% y = [u v eta w]; u = r*lambda2; v = r*lambda2'; w = eta'; Ref: Patil and
% DasGupta 2013

%     tic
    options = odeset('stats','off','Mass',@mass,'RelTol',1e-5,'AbsTol',1e-13);
    [rvals, sol] = ode45(@fnIVP,rrange,y0,options);
%     toc
%     r
%     y

%     [rvals, sol] = RK4manual(h, P, MU, mus, alphas, y0, rrange);


    n = size(sol,1);
    m = size(rvals,1);
    if n~=m
       disp('Check. lengths of r and y are not same here!'); 
    end
    rhoat1 = sol(n,1);
    etaat1 = sol(n,3);
    

    function dydt = fnIVP(r,y)
        
        if r == 0
           lambda2 = lambda0; 
        else
            lambda2 = y(1)/r;
        end
        
        [~, ~, CC, ~, ~, FF] = ode45essentials(y.', lambda2, r, h, P, chixMU, mus, alphas);
        dydt = [lambda2 + y(2); CC; y(4); FF];
        
    end

    function M = mass(r,y)
        
        if r == 0
           lambda2 = lambda0; 
        else
            lambda2 = y(1)/r;
        end
        
        [AA, BB, ~, DD, EE, ~] = ode45essentials(y.', lambda2, r, h, P, chixMU, mus, alphas);
        M = [1 0 0 0; 0 AA 0 BB; 0 0 1 0; 0 DD 0 EE];
        
    end

end