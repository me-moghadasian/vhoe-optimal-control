function res = doInterp(M, n, PROBLEM)
    colN = size(M, 2);
    tN   = numel(PROBLEM.time);
    res  = zeros(n*PROBLEM.c, colN);
    for j = 1:colN
        for k = 1:n
            if size(M, 1)==n % COnstant time matrix
                res((k-1)*PROBLEM.c+1:k*PROBLEM.c, j) = M(k, j);
            else
                res((k-1)*PROBLEM.c+1:k*PROBLEM.c, j) =...
                    evalAtAdjNodes(M((k-1)*tN+1:k*tN, j),...
                                   PROBLEM.time, PROBLEM.adjTime);
            end
        end
    end
end