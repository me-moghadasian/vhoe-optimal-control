function zR = boundaryReverse(boundary, zA)
    a = boundary(1);
    b = boundary(2);
    zR = 2*(zA-a)./(b-a)-1;    
end