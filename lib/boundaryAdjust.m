function zA = boundaryAdjust(boundary, z)
    a = boundary(1);
    b = boundary(2);
    zA = (z+1).*(b-a)./2+a;
end

