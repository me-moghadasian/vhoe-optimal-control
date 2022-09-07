function delu = calcDelU(A, nx, m, noT)
n = 2*nx;
W = calcW(n, m);
% W = zeros(n, rCombRep(n, m));
% 
% i1 = 1;
% p  = m;
% q  = p;
% i2 = 1;
% 
% 
% while i1~=(n+1)
%     if ((i1~=1)&&(p==-1))||(i1==n) 
%         p = m-sum(W(1:i1-1, i2));
%         q = p;
%     end 
%     r  = n-i1;
%     c = rCombRep(r, q-p);
%     W(i1,i2:i2+c-1) = p;
%     i2 = i2+c;
%     p = p-1;
%     if i2==(rCombRep(n, m)+1)
%         i1 = i1+1;
%         i2 = 1;
%     end
% end

delu = zeros(nx*noT, rCombRep(n, m)*noT);
c = 1;
for k = 1:rCombRep(n, m)
    sind = (k-1)*noT+1;
    find =  k*noT;
    if isz(W(nx+1:end, k))
        delu(:, sind:find) = makeItDiag(A(:, c), nx, noT);
        c = c + 1;
    end
end
end
function r = isz(v)
    r = 1;
    for k = 1:numel(v)
        if v(k)~=0
            r = 0;
            break;
        end
    end
end
