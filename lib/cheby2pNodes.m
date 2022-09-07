function z = cheby2pNodes(m)
    z = zeros(m, 1);
    for k = 1:m
        z(k) = cos((k-1)*pi/(m-1)); %(n-1)th second kind plus +-1
    end
    
end

