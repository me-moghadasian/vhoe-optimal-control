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
PROBLEM.type = types{1};

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
C.ic = ic;
C.Q = Q;
C.R = R;
C.B = B;

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
C.q = 1;
C.Gain = [];
C.Gain = PROBLEM.S{1}(PROBLEM.c*PROBLEM.n+1:PROBLEM.c:end,:);

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true), tspan, ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 1st order: ' num2str(J(end))])

plot(x, v, 'DisplayName','1st order','LineWidth',1,'LineStyle','-');

% with controller 3nd order
C.q = 3;
C.Gain = [];
for k=1:C.q
    C.Gain = [C.Gain,...
        PROBLEM.S{k}(PROBLEM.c*PROBLEM.n+1:PROBLEM.c:end,:)];
end

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true), tspan, ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 3rd order: ' num2str(J(end))])

plot(x, v, 'DisplayName','3rd order','LineWidth',1,'LineStyle','-');
legend show

% with controller 5nd order
C.q = 5;
C.Gain = [];
for k=1:C.q
    C.Gain = [C.Gain,...
        PROBLEM.S{k}(PROBLEM.c*PROBLEM.n+1:PROBLEM.c:end,:)];
end

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true), tspan, ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 5th order: ' num2str(J(end))])

plot(x, v, 'DisplayName','5th order','LineWidth',1,'LineStyle','-');
legend show

% with controller 5nd order
C.q = 7;
C.Gain = [];
for k=1:C.q
    C.Gain = [C.Gain,...
        PROBLEM.S{k}(PROBLEM.c*PROBLEM.n+1:PROBLEM.c:end,:)];
end

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[controlOnTime, controlOnStates] =...
    ode45(@(t, states)dyEqu(t, states, true), tspan, ic, opts);

x = controlOnStates(:, 1);
v = controlOnStates(:, 2);
J = controlOnStates(:, 3);
disp(['performance index for 7th order: ' num2str(J(end))])

plot(x, v, 'DisplayName','7th order','LineWidth',1,'LineStyle','-');
legend show

saveas(f, 'traj.fig')
