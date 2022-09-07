function res = rCombRep(n, r)
    if n==0
        res = 1;
    else
        res = factorial(n+r-1)/factorial(r)/factorial(n-1);
    end
end