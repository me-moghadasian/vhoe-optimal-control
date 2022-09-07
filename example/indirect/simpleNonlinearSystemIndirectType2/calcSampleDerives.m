function [t, A, B, Q, R, n, c, mp, q] = calcSampleDerives()
    global C
    n = 2;
    mp = 1;
    c = 15;
    q = 9;

    t = [0 1]';
    ov = ones(size(t));

    alp = C.alp;
    bet = C.bet;
    del = C.del;
    
    %dx/dt = v;
    %dv/dt = alp*x+del*v+bet*x^3+u;
    
    A{1} = [0*ov,     1*ov;
            alp*ov, del*ov];
    A{2} = [0*ov, 0*ov, 0*ov;
            0*ov, 0*ov, 0*ov];
    A{3} = [  0*ov, 0*ov, 0*ov, 0*ov;
            bet*ov, 0*ov, 0*ov, 0*ov]; 
    A{4} = zeros(2*numel(t), rCombRep(2, 4));
    A{5} = zeros(2*numel(t), rCombRep(2, 5));
    A{6} = zeros(2*numel(t), rCombRep(2, 6));
    A{7} = zeros(2*numel(t), rCombRep(2, 7));
    A{8} = zeros(2*numel(t), rCombRep(2, 8));
    A{9} = zeros(2*numel(t), rCombRep(2, 9));
    
    B = [0*ov;
         1*ov];  
        
    Q = [1*ov, 0*ov;
         0*ov, 1*ov];
     
    R = [1*ov]; 
end
 
 
 






