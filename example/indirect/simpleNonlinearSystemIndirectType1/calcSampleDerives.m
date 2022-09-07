function [t, A, B, Q, R, n, c, mp, q] = calcSampleDerives()
    global C
    n = 2;
    mp = 1;
    q = 7;

    t = [0:0.25:1]';
    c = numel(t);
    
    alp = C.alp;
    bet = C.bet;
    del = C.del;
    
    %dx/dt = v;
    %dv/dt = alp*x+del*v+bet*x^3+u;
    
    A{1} = [0,     1;
            alp, del];
    A{2} = [0, 0, 0;
            0, 0, 0];
    A{3} = [  0, 0, 0, 0;
            bet, 0, 0, 0]; 
    A{4} = zeros(2, rCombRep(2, 4));
    A{5} = zeros(2, rCombRep(2, 5));
    A{6} = zeros(2, rCombRep(2, 6));
    A{7} = zeros(2, rCombRep(2, 7));
    A{8} = zeros(2, rCombRep(2, 8));
    A{9} = zeros(2, rCombRep(2, 9));
    
    B = [0;
         1];  
        
    Q = [1, 0;
            0, 1];
     
    R = 0.05*[1]; 

end
 
 
 






