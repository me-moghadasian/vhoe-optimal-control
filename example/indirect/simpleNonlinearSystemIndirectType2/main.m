clc
clear all

tic
global C

%% Problem specific paramters initialization
C.alp = -1;
C.bet = -5;
C.del = -0.02;

%% Read problem data
types{1} = 'InfiniteHorizon';
types{2} = 'FiniteHorizon-FixedTerminalTime-ZeroTerminalState';
types{3} = 'FiniteHorizon-FixedTerminalTime-MayerTerminalState';

PROBLEM.mode = 'indirect'; 
PROBLEM.type = types{2};

[t, A, B, Q, R, n, c, m, q] = calcSampleDerives();

PROBLEM.n    = n;
PROBLEM.time = t;
PROBLEM.c    = c;
PROBLEM.m    = m;
PROBLEM.q    = q;

PROBLEM.A = A;
PROBLEM.B = B;
PROBLEM.Q = Q;
PROBLEM.R = R;
PROBLEM.H = zeros(n);
%
PROBLEM.D = [];

%% Solve problem
PROBLEM = VHOE_NOCSolver(PROBLEM);

save Result PROBLEM
toc

%% Testing the solution
tspan = [0 100];
ic = [1; 0; 0];
C.Q = Q([1:numel(PROBLEM.time):end],:);
C.R = R([1:numel(PROBLEM.time):end],:);

% without controller
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[controlOffTime, controlOffStates] =...
    ode45(@(t, states)dyEqu(t, states, false), tspan, ic, opts);

x = controlOffStates(:, 1);
v = controlOffStates(:, 2);

f = figure();
plot(x, v, 'DisplayName','without Control','LineWidth',1,'LineStyle','-');
xlabel('x')
ylabel('v')
hold on

% with controller 1rd order
perVec = makePerVec(ic(1:2), 1);

S_aug = [PROBLEM.S{1}];
Z = S_aug*perVec;
p1_aug = Z((PROBLEM.n+0)*PROBLEM.c+1:(PROBLEM.n+0)*PROBLEM.c+PROBLEM.c);
p2_aug = Z((PROBLEM.n+1)*PROBLEM.c+1:(PROBLEM.n+1)*PROBLEM.c+PROBLEM.c);

C.u = -p2_aug./PROBLEM.R;
C.t = PROBLEM.STime;

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true),...
    [PROBLEM.STime(1), PROBLEM.STime(end)], ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 1st order: ' num2str(J(end))])

plot(x, v, 'DisplayName','1st order','LineWidth',1,'LineStyle','-');

% with controller 3rd order
perVec = makePerVec(ic(1:2), 3);

S_aug = [PROBLEM.S{1}, PROBLEM.S{2}, PROBLEM.S{3}];
Z = S_aug*perVec;
p1_aug = Z((PROBLEM.n+0)*PROBLEM.c+1:(PROBLEM.n+0)*PROBLEM.c+PROBLEM.c);
p2_aug = Z((PROBLEM.n+1)*PROBLEM.c+1:(PROBLEM.n+1)*PROBLEM.c+PROBLEM.c);

C.u = -p2_aug./PROBLEM.R;
C.t = PROBLEM.STime;

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true),...
    [PROBLEM.STime(1), PROBLEM.STime(end)], ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 3rd order: ' num2str(J(end))])

plot(x, v, 'DisplayName','3rd order','LineWidth',1,'LineStyle','-');

% with controller 5th order
perVec = makePerVec(ic(1:2), 5);

S_aug = [PROBLEM.S{1}, PROBLEM.S{2}, PROBLEM.S{3},...
    PROBLEM.S{4}, PROBLEM.S{5}];
Z = S_aug*perVec;
p1_aug = Z((PROBLEM.n+0)*PROBLEM.c+1:(PROBLEM.n+0)*PROBLEM.c+PROBLEM.c);
p2_aug = Z((PROBLEM.n+1)*PROBLEM.c+1:(PROBLEM.n+1)*PROBLEM.c+PROBLEM.c);

C.u = -p2_aug./PROBLEM.R;
C.t = PROBLEM.STime;

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true),...
    [PROBLEM.STime(1), PROBLEM.STime(end)], ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 5th order: ' num2str(J(end))])

plot(x, v, 'DisplayName','5th order','LineWidth',1,'LineStyle','-');

% with controller 7th order
perVec = makePerVec(ic(1:2), 7);

S_aug = [PROBLEM.S{1}, PROBLEM.S{2}, PROBLEM.S{3},...
    PROBLEM.S{4}, PROBLEM.S{5}, PROBLEM.S{6}, PROBLEM.S{7}];
Z = S_aug*perVec;
p1_aug = Z((PROBLEM.n+0)*PROBLEM.c+1:(PROBLEM.n+0)*PROBLEM.c+PROBLEM.c);
p2_aug = Z((PROBLEM.n+1)*PROBLEM.c+1:(PROBLEM.n+1)*PROBLEM.c+PROBLEM.c);

C.u = -p2_aug./PROBLEM.R;
C.t = PROBLEM.STime;

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true),...
    [PROBLEM.STime(1), PROBLEM.STime(end)], ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 7th order: ' num2str(J(end))])

plot(x, v, 'DisplayName','7th order','LineWidth',1,'LineStyle','-');
legend show

saveas(f, 'traj.fig')
