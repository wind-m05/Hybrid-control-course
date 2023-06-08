clear all
close all
clc
%%
% Initial condition
x0 = [1; 0];

% Simulation horizon
% The simulation stops when the end of either TSPAN or JSPAN has been reached
TSPAN = [0 10]; % simulation time interval
JSPAN = [0 20]; % simulation interval for discrete jumps

% Rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 2;

% ode45 solver options
RelTol = 1e-6;
MaxStep = 0.1;
options = odeset('RelTol',RelTol,'MaxStep',MaxStep);

% Simulate using Hybrid Equations Toolbox
[t,j,x] = HyEQsolver(@a_f, @a_g, @a_C, @a_D, x0,TSPAN,JSPAN,rule,options);


% Plot solution over time
figure
subplot(2,1,1)
plotflows(t,j,x(:,1))
grid on
xlabel('time t')
ylabel('x_1(t)')

subplot(2,1,2)
plotflows(t,j,x(:,2))
grid on
xlabel('time t')
ylabel('x_2(t)')
