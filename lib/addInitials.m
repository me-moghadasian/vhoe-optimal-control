function res = addInitials(sol, iniTotal, noVar, m)
    res = zeros(m, noVar);
    for i1 = 1:noVar
        res(:, i1) = [sol([1:m-1]+(i1-1)*(m-1), 1); iniTotal(i1)];
    end
end