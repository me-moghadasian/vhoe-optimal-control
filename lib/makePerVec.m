function perVec = makePerVec(v, m)
    n      = numel(v);
    perVec = [];
    for k = 1:m
        W = calcW(n, k);
        pv = ones(rCombRep(n, k), 1);
        for j = 1:rCombRep(n, k)
            for i = 1:n
                pv(j, 1) = pv(j, 1)*v(i)^W(i, j)/factorial(W(i, j));
%                 pv(j, 1) = pv(j, 1)*v(i)^W(i, j);
            end
        end
        pv = pv*factorial(k);
        perVec = [perVec; pv];
    end
end

function W = calcW(n, m)
    W = zeros(n, rCombRep(n, m));

    i1 = 1;
    p  = m;
    q  = p;
    i2 = 1;

    while i1~=(n+1)
        if ((i1~=1)&&(p==-1))||(i1==n) 
            p = m-sum(W(1:i1-1, i2));
            q = p;
        end 
        r  = n-i1;
        c = rCombRep(r, q-p);
        W(i1,i2:i2+c-1) = p;
        i2 = i2+c;
        p = p-1;
        if i2==(rCombRep(n, m)+1)
            i1 = i1+1;
            i2 = 1;
        end
    end
end

function res = rCombRep(n, r)
    if n==0
        res = 1;
    else
        res = factorial(n+r-1)/factorial(r)/factorial(n-1);
    end
end