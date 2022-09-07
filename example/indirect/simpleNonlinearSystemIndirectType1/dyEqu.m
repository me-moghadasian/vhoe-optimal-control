function dstates = dyEqu(t, states, isControlOn)
    global C
    alp = C.alp;
    bet = C.bet;
    del = C.del;
    
    x = states(1);
    v = states(2);
    
    if isControlOn
        s = [x; v];
        pVec = makePerVec(s, C.q);
        pVec(1:2) = s;
        u = -C.R^-1*C.B'*C.Gain*pVec;
        dj = 1/2*s'*C.Q*s+1/2*u'*C.R*u;
    else 
        u = 0;
        dj = 0;
    end
    
    
    dxdt = v;
    dvdt = alp*x+del*v+bet*x^3+u;
    
    dstates = [dxdt; dvdt; dj];
end