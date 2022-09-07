function [sensis, opt] = calcSensis(PROBLEM, k, opt, C)
    eachOrdCount = rCombRep(PROBLEM.n, k);
    sensis       = zeros(PROBLEM.c*PROBLEM.n*2, eachOrdCount);
    fnl = zeros(PROBLEM.n, 1); 
    ini = zeros(PROBLEM.n, 1);
    
    for j = 1:eachOrdCount
        ini = ini*0;
        if k==1
            ini(j) = 1;
        end
        
        if strcmp(PROBLEM.type, C.types{1})
            [~, sol, opt] = linIHSolver(PROBLEM.A{1}, PROBLEM.B, PROBLEM.Q,...
                PROBLEM.R, PROBLEM.nonH{k}(:, j),...
                PROBLEM.time, ini, opt);
        elseif strcmp(PROBLEM.type, C.types{2})
            [~, sol, opt] = linFixedTerminalStateBVSolver(PROBLEM.D{1, 1},...
                PROBLEM.nonH{k}(:, j), PROBLEM.ts, ini, fnl, opt);
        elseif strcmp(PROBLEM.type, C.types{3})
            [~, sol, opt] = linMayerTerminalStateBVSolver(PROBLEM.D{1, 1},...
                PROBLEM.nonH{k}(:, j), PROBLEM.H, PROBLEM.ts, ini, opt);
        else
            disp('Type error');
            return;
        end
        
        sensis(:, j) = reshape(sol, numel(sol), 1);
    end
end