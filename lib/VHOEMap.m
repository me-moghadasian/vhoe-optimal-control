function PHI = VHOEMap(n, I, J)
%% Initialization
PHI = cell(I, J);
PHI{1, 1} = (1:n)';

%% Mappings
for i_ = 2:I
    for j_ = 1:J
        if j_==1
            %% Simple mapping
            PHI = SMapping(n, i_, j_, PHI);
        elseif (i_==j_+1)&&(j_~=1)
            %% Full mapping            
            PHI = FMapping(n, i_, j_, PHI);
        elseif (i_>j_+1)&&(j_~=1)
            %% Partial mapping
            PHI = PMapping(n, i_, j_, PHI);
        end
    end
end
end
function PHI = SMapping(n, i_, j_, PHI)
%% Simple mapping
ir = 1;
ic = 1;
PHI{i_, j_} = zeros(rCombRep(n,i_), rCombRep(n, i_-j_));
for k=1:n
    rp  = xi(i_-1, k, n);
    cp  = xi(i_-2, k, n);
    irp = rCombRep(n, i_-1) - rp + 1;
    icp = rCombRep(n, i_-2) - cp + 1;
    phi = PHI{i_-1, j_}(irp:end,icp:end);
    PHI{i_, j_}(ir:ir+rp-1,ic:ic+cp-1) = phi;
    for p = 1:(size(PHI{i_, j_}, 2)-ic)
        PHI{i_, j_}(ir+p,ic+p) = phi(1, 1);
    end
    ir = ir + rp;
    ic = ic + cp;
end
end
function PHI = FMapping(n, i_, j_, PHI)
%% Fully compressed mapping
PHI{i_, j_} = zeros(rCombRep(n,i_), rCombRep(n, i_-j_));

ir = 1;
ic = 1;
v  = 0;
for k=1:n
    rp  = xi(i_-1, k, n);
    cp  = xi(1, k, n);
    irp = rCombRep(n, i_-1) - rp + 1;
    icp = rCombRep(n, 1)    - cp + 1;
    phi = PHI{j_, j_-1}(irp:end,icp:end);
    
    exind = xi(j_-1, k, n);
    for p = exind+1:size(phi, 1)
        phi(p, 1) = phi(p-1, 1) + 1;
    end
    phi = phi + v.*(phi~=0);
    PHI{i_, j_}(ir:ir+rp-1,ic:ic+cp-1) = phi;
    
    ir = ir + rp;
    ic = ic + 1;
    if k<n
        v  = v + phi(exind+1, 1) - phi(exind+1, 2);
    end
end
end
function PHI = PMapping(n, i_, j_, PHI)
%% Partial mapping

ir = 1;
ic = 1;
PHI{i_, j_} = zeros(rCombRep(n,i_), rCombRep(n, i_-j_));
for k=1:n
    rp  = xi(i_-1,    k, n);
    cp  = xi(i_-j_-1, k, n);
    irp = rCombRep(n, i_-1)    - rp + 1;
    icp = rCombRep(n, i_-j_-1) - cp + 1;
    phi = PHI{i_-1, j_}(irp:end,icp:end);
    PHI{i_, j_}(ir:ir+rp-1,ic:ic+cp-1) = phi;
    for p = 1:(size(PHI{i_, j_}, 2)-ic)
        PHI{i_, j_}(ir+p,ic+p) = phi(1, 1);
    end
    
    PHI{i_, j_} = recurTreeGen(n, i_, j_, j_-1,...
        PHI, k, k+1, PHI{i_, j_}, ir, ic);
    
    ir = ir + rp;
    ic = ic + cp;
end
end
function res = recurTreeGen(n, i_, j_, bNo, PHI, lev, depth, mem, ir, ic)
    % initialization
    if (lev >= n)||(depth > n)
        res = mem;
        return;
    end
    
    % all branch in a level loop
    for k=1:bNo
        irp_f = rCombRep(n, i_-j_-1+k);
        icp_f = rCombRep(n, i_-j_-1);
        rp    = xi(i_-j_-1+k, depth, n);
        cp    = xi(i_-j_-1,   depth, n);
        irp_i = irp_f - rp + 1;
        icp_i = icp_f - cp + 1;
        phi   = PHI{i_-1, j_}(irp_i:irp_f,icp_i:icp_f);       
        irn = ir + xi(i_-j_-1+k, depth-1, n);
        icn = ic + xi(i_-j_-1,   depth-1, n);
        
        v = findAugValue(mem, phi(1, 1), irn, icn);
        phi = phi + v.*(phi~=0);
        
        mem(irn:irn+rp-1, icn:icn+cp-1) = phi;
        for p = 1:(size(mem, 2)-icn)
            mem(irn+p,icn+p) = phi(1, 1);
        end
        
        % start new depth
        if (depth < n)&&(lev <= n)
            mem = recurTreeGen(n, i_, j_, k, PHI, lev, depth+1,...
                mem, irn, icn);
        end
    end
    
    % return data
    res = mem;
end
function v   = findAugValue(mem, phiH, irn, icn)
    v = 0;
    for k= irn-1:-1:1
        if mem(k, icn)~=0
            v = mem(k, icn) - phiH + 1;
            break;
        end
    end
end





