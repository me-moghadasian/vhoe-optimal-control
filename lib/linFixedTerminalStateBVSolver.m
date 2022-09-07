function [time, sol, opt] = linFixedTerminalStateBVSolver(A, nonH, boundary, ini, fnl, opt)
    ini = makeItClmn(ini);
    fnl = makeItClmn(fnl);
    noVar = numel(ini)*2;
    index_ini = zeros(1, noVar); 
    m = size(A, 1)/noVar;
    a = boundary(1);
    b = boundary(2);
    D = cheby2pD(m);
    Dprime = 2/(b-a)*D;
    augD = [];
    for k = 1:noVar
        augD = blkdiag(augD, Dprime);
        index_ini(k) = m+(k-1)*m; 
    end
    G = augD-A;
    opt.G = G;
    F = G;
    F(index_ini, :) = [];
    nonH(index_ini, :) = [];
    F_ini = F(:, index_ini);
    F(:, index_ini) = [];
    if isempty(opt.FInv)
        Phi = F^-1;
        opt.FInv = Phi;
    else
        Phi = opt.FInv;
    end
    Phi_UL = Phi(1:end/2, 1:end/2);
    Phi_UR = Phi(1:end/2, end/2+1:end);
    nonH_U = nonH(1:end/2, 1);
    nonH_L = nonH(end/2+1:end, 1);
    F_ini_UL = F_ini(1:end/2, 1:end/2);
    F_ini_LR = F_ini(end/2+1:end, end/2+1:end);
    F_ini_UR = F_ini(1:end/2, end/2+1:end);
    F_ini_LL = F_ini(end/2+1:end, 1:end/2);
    M_known = (Phi_UL(1:m-1:end, :)*F_ini_UL+...
        Phi_UR(1:m-1:end, :)*F_ini_LL)*ini;
    Mp_unknown = Phi_UR(1:m-1:end, :)*F_ini_LR+...
        Phi_UL(1:m-1:end, :)*F_ini_UR;
    nonHp_U = Phi_UL(1:m-1:end, :)*nonH_U;
    nonHp_L = Phi_UR(1:m-1:end, :)*nonH_L;
    Hp_total = nonHp_U + nonHp_L;
    costate_ini = -Mp_unknown\(fnl+M_known-Hp_total);
    iniTotal = [ini; costate_ini];
    sol = -Phi*F_ini*iniTotal+Phi*nonH;
    sol = addInitials(sol, iniTotal, noVar, m);
    cNodes = cheby2pNodes(m);
    time = boundaryAdjust(boundary, cNodes);
end 














