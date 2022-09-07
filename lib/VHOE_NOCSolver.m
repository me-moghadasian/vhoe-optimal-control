function PROBLEM = VHOE_NOCSolver(PROBLEM)
    %% VHOE_NOCSolver is designed to solve nonlinear optimal
    % control problems. This program uses Vectorised High
    % Order Expansions (VHOE) method which was introduced
    % by Dr. Mehdi Moghadasian in 2019. 
    %
    % This tool supports 3 types of optimal control problems:
    %    1. Infinite Horizon 
    %    2. Finite Horizon with Fixed Terminal Time and State
    %       (Hard Constraint Problem)
    %    3. Finite Horizon with Fixed Terminal Time and Mayer
    %       Termial State (Soft Constraint Problem)
    %
    % Check \examples\ folder or the following reference 
    %  for more information.
    % [1] ...
    %
    % -------------------------------------------------------
    % The variable PROBLEM contains both input and output of
    %   the function.
    %
    % Inputs:
    %	PROBLEM.n: Number of states
    %   PROBLEM.time: A vector that contains all time instances at
    %       which the coefficients have been evaluated by user for
    %       types 2 and 3 problems;
    %       or the time vector that should be used for
    %       forward and backward Euler integration fo type 1 problem.
    %   PROBLEM.c: Number of collocation points for pseudo-spectral method
    %       which is used to solve TPBV problem for type 2 and 3; or       
    %   PROBLEM.m: Number of control variables.
    %   PROBLEM.q: Solution order
    %   PROBLEM.A: Mathematical model coefficients related to states
    %   PROBLEM.B: Mathematical model coefficients related to control
    %       variables
    %   PROBLEM.Q: Weight matrix in quadratic format related to state
    %       variables
    %   PROBLEM.R: Weight matrix in quadratic format related to control
    %       variables
    %   PROBLEM.H: Weight matrix in quadratic format related to state
    %       variables in Mayer form
    %   PROBLEM.D: Combined state-costate coefficients if the provided
    %       model is prepared explicitly by applying OC theory.
    %   
    %   *Note: For types 2 and 3, all coefficients are time dependent; so
    %       each elements of the matrices must be provided as a vactor 
    %       even if the coeffs are constant. For type 1, all elements are
    %       constant scalar.
    %
    % Outputs:
    %   
    
    %% Check problem data
    % The problem type will be specified here
    C.types{1} = 'InfiniteHorizon';
    C.types{2} = 'FiniteHorizon-FixedTerminalTime-ZeroTerminalState';
    C.types{3} = 'FiniteHorizon-FixedTerminalTime-MayerTerminalState';
    
    % The values t0, tf and ts are assined here in struct
    %   PROBLEM. This struct is entended to store program specific data
    % t0: inital time
    % tf: final time
    % ts: time range
    if strcmp(PROBLEM.type, C.types{1})
        % For infinite horizon problem t0 and tf are not defined
        PROBLEM.t0 = 0;
        PROBLEM.tf = inf;
    else
        PROBLEM.t0 = PROBLEM.time(1);
        PROBLEM.tf = PROBLEM.time(end);
    end
    PROBLEM.ts = [PROBLEM.t0, PROBLEM.tf];

    %% Interpolation on collocation points
    PROBLEM.adjTime = [];
    if ~strcmp(PROBLEM.type, C.types{1})
        PROBLEM.cNodes  =  cheby2pNodes(PROBLEM.c);
        PROBLEM.adjTime = boundaryAdjust(PROBLEM.ts, PROBLEM.cNodes);
        if strcmp(PROBLEM.mode, 'indirect')
            for k = 1:PROBLEM.q
                PROBLEM.A{k} = doInterp(PROBLEM.A{k}, PROBLEM.n, PROBLEM);
            end
            PROBLEM.B = doInterp(PROBLEM.B, PROBLEM.n, PROBLEM);
            PROBLEM.Q = doInterp(PROBLEM.Q, PROBLEM.n, PROBLEM);
            PROBLEM.R = doInterp(PROBLEM.R, PROBLEM.m, PROBLEM);
        else
            for k = 1:PROBLEM.q
                PROBLEM.D{k} = doInterp(PROBLEM.D{k}, 2*PROBLEM.n, PROBLEM);
            end
        end
    end
    if strcmp(PROBLEM.type, C.types{1})
        PROBLEM.c = numel(PROBLEM.time);
    end
    
    %% Solve problem
    PROBLEM.PHIx = VHOEMap(1*PROBLEM.n, PROBLEM.q, PROBLEM.q-1);
    PROBLEM.PHIz = VHOEMap(2*PROBLEM.n, PROBLEM.q, PROBLEM.q-1);
    PROBLEM.D    = calcDelta(PROBLEM, C.types);

    PROBLEM.S    = cell(PROBLEM.q, 1);
    PROBLEM.nonH = cell(PROBLEM.q, 1);
    PROBLEM.Gam  = cell(PROBLEM.q, PROBLEM.q);
    
    opt.FInv      = [];
    opt.KBar      = [];
    opt.TCoeff = [];

    for k = 1:PROBLEM.q
        % Calc non-homogenous part
        [PROBLEM.nonH{k}, Gam] = calcNonHomo(PROBLEM, k);
    
        % Calc sensitivities
        [PROBLEM.S{k}, opt] = calcSensis(PROBLEM, k, opt, C);
    
        % Augmenting Gamma
        PROBLEM.Gam = augmentGamma(PROBLEM, k, Gam);
    end

    %% Post-processing
    if ~strcmp(PROBLEM.type, C.types{1})
        for k = 1:PROBLEM.q
            for j = 1:2*PROBLEM.n
                sind = (j-1)*PROBLEM.c+1;
                find = j*PROBLEM.c;
                PROBLEM.S{k}(sind:find,:) =...
                    flipud(PROBLEM.S{k}(sind:find,:));
            end
        end
        PROBLEM.STime = flipud(PROBLEM.adjTime);
    else
        PROBLEM.STime = PROBLEM.time;
    end
end