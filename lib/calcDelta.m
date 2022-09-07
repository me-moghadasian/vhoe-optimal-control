function D = calcDelta(PROBLEM, types)
    D = cell(PROBLEM.q, 1);
    if strcmp(PROBLEM.mode,'direct')
        for k = 1:PROBLEM.q
            D{k, 1}  = makeItDiag(PROBLEM.D{k, 1},...
                2*PROBLEM.n, PROBLEM.c);
        end
        return;
    end
    
    A = cell(PROBLEM.q, 1);
    if ~strcmp(PROBLEM.type, types{1})
        A    = PROBLEM.A;
        B    = PROBLEM.B;
        Q    = PROBLEM.Q;
        R    = PROBLEM.R;
    else
        for qc = 1:PROBLEM.q
            A{qc} = repelem(PROBLEM.A{qc}, PROBLEM.c, 1);
        end
        B    = repelem(PROBLEM.B, PROBLEM.c, 1);
        Q    = repelem(PROBLEM.Q, PROBLEM.c, 1);
        R    = repelem(PROBLEM.R, PROBLEM.c, 1);
    end
    n    = PROBLEM.n;
    PHIx = PROBLEM.PHIx;
    
    A1d = makeItDiag(A{1}, PROBLEM.n, PROBLEM.c);
    Bd  = makeItDiag(B,    PROBLEM.n, PROBLEM.c);
    Qd  = makeItDiag(Q,    PROBLEM.n, PROBLEM.c);
    Rd  = makeItDiag(R,    PROBLEM.m, PROBLEM.c);
    D{1, 1} = [ A1d, -Bd*Rd^-1*Bd'
               -Qd,    -A1d'];
    for k = 2:PROBLEM.q
        delu  = calcDelU(A{k}, n, k, PROBLEM.c);
        dell  = calcDelL(PHIx{k, 1}, -A{k}, n, k, PROBLEM.c);
        D{k, 1}  = [delu; dell];
    end
end