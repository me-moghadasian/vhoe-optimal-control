function r = xi(i, k, n)
    N = zeros(n, 1);
    N(n, 1) = 1;
    for j = 1:i
        for p = 1:n
            N(p) = sum(N(p:end));
        end
    end
    r = sum(N(k:end));
end