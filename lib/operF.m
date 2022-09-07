function res = operF(M, PHI, I, J, noT)
    res = [];
    for j_=1:J 
        sind = (j_-1)*noT+1;
        find =    j_*noT;
        res = [res putInFormat(M(:, sind:find), PHI{I, 1}, noT)];
    end
end