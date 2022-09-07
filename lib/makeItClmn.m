function res = makeItClmn(v)
    res = v;
    if size(v, 2)~=1
        res = v';
    end
end