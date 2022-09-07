function w = evalAtAdjNodes(F, x, xA)
    m = size(xA, 1);
    w = zeros(m, 1);
    for i1 = 1:m
        %%w(i1, 1) = interpn(x, F, xA(i1),'spline');
		w(i1, 1) = interpn(x, F, xA(i1));
    end
end

