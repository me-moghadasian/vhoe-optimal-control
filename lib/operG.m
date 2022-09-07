function res = operG(M, PHI, I, J, noT)
    noRowData = size(M, 1)/noT;
    r = size(PHI{I, J}, 1);
    c = size(PHI{I, J}, 2);
    res = sparse(c*noT*noRowData, r*noT); % c and r are replaced because transpose is needed
    for i1 = 1:r
        for i2 = 1:c
            ind = PHI{I, J}(i1, i2);
            if (ind ~= 0)
                res((i2-1)*noT*noRowData+1:i2*noT*noRowData,...
                    (i1-1)*noT+1:i1*noT) =...
                    M(:,(ind-1)*noT+1:ind*noT);
            end
        end
    end
end