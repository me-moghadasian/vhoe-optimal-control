function [time, sol, opt] = linIHSolver(A, B, Q, R, nonH, time, ini, opt)
    c = numel(time);
    n = numel(ini);
    
    if isempty(opt.KBar)
        [KBar,~,~] = care(A,B,Q,R);
        opt.KBar = KBar;
    else
        KBar = opt.KBar;
    end
    
    if isempty(opt.TCoeff)
        TCoeff     = KBar*(-B*R^-1*B')-(-A');
        opt.TCoeff = TCoeff;
    else
        TCoeff = opt.TCoeff;
    end
    
    nonHx = reshape(nonH(1:c*n), c, n);
    nonHl = reshape(nonH(c*n+1:end), c, n);
    
    sol = zeros(c, 2*n);
    
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [~, T] = ode45(@(t, y) dequTm(t, y, time, TCoeff, KBar, ...
        full(nonHx), full(nonHl)),...
        time, zeros(1, n), opts);
    T = flipud(T);
    
    [~, sx] = ode45(@(t, y) dequsx(t, y, time, T, KBar, full(nonHx),...
        A, B, R),...
        time, ini, opts);
    sol(:, 1:n) = sx;
    
    for ct = 1:c
        sol(ct, n+1:end) = (KBar*sol(ct, 1:n)'+T(ct, :)')';
    end
    
    sol = reshape(sol, numel(sol), 1);
end

function dTm = dequTm(tm, T, time, TCoeff, KBar, nonHx, nonHl)
    tc = time(end) - tm;
    nonHlc = interp1(time, nonHl, tc);
    nonHxc = interp1(time, nonHx, tc);
    dTm = -(-TCoeff*T+...
            (nonHlc'-KBar*nonHxc'));
end

function dsx = dequsx(t, sx, time, T, KBar, nonHx, A, B, R)
    nonHxc = interp1(time, nonHx, t);
    Tc = interp1(time, T, t);
    sl = KBar*sx+Tc';
    u = -R^-1*B'*sl;
    dsx =  A*sx + B*u + nonHxc';
end
