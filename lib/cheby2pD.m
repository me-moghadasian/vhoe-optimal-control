function D = cheby2pD(m)
    D = zeros(m);
    y = cheby2pNodes(m);
    n = m-1;
    for j = 0:n
        for k = 0:n
            if ((j==0)&&(k==0))
                D(j+1, k+1) = (1+2*n^2)/6;
            elseif ((j==n)&&(k==n))
                D(j+1, k+1) = -(1+2*n^2)/6;
            elseif ((j==0)&&(k==n))
                D(j+1, k+1) = 1/2*(-1)^n;
            elseif ((j==n)&&(k==0))
                D(j+1, k+1) = -1/2*(-1)^n;
            elseif ((j==0)&&(k>0)&&(k<n))
                D(j+1, k+1) = 2*(-1)^k/(1-y(k+1));
            elseif ((k==0)&&(j>0)&&(j<n))
                D(j+1, k+1) = -1/2*(-1)^j/(1-y(j+1));
            elseif ((j==n)&&(k>0)&&(k<n))
                D(j+1, k+1) = -2*(-1)^(n-k)/(1+y(k+1));
            elseif ((k==n)&&(j>0)&&(j<n))
                D(j+1, k+1) = 1/2*(-1)^(n-j)/(1+y(j+1));
            elseif ((j==k)&&(j>0)&&(j<n)) 
                D(j+1, k+1) = -1/2*y(j+1)/(1-y(j+1)^2);
            else
                D(j+1, k+1) = (-1)^(k-j)/(y(j+1)-y(k+1));
            end
        end
    end
end