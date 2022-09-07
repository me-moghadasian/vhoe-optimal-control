function res = putInFormat(data, PHI, noT)
    if isempty(PHI)
        res = [];
        return
    end
    
    r = size(PHI, 1);
    c = size(PHI, 2);
    b = size(data, 2)/noT; % batch size
    res = sparse(r*noT, c*noT*b);
    for i1 = 1:r
        for i2 = 1:c
            ind = PHI(i1, i2);
            if (ind ~= 0)
                for i3 = 1:b
                    sri = (i1-1)*noT+1;
                    fri = sri-1+noT;
                    sci = (i2-1)*b*noT+(i3-1)*noT+1;
                    fci = sci-1+noT;
                        
                    res(sri:fri, sci:fci) = ...
                            data((ind-1)*noT+1:ind*noT,...
                                  (i3-1)*noT+1:i3*noT);
                        
%                     res(sri:fri, sci:fci) = ...
%                             diag(data((ind-1)*noT+1:ind*noT, i3));
                end
%                     res((i1-1)*noT+1:i1*noT,...
%                         (i2-1)*noT+1:i2*noT) = ...
%                             diag(data((ind-1)*noT+1:ind*noT, :));
            end
        end
    end
end

