clear
close all
clc

set(0,'defaulttextInterpreter','latex')

%% Initial conditions and timespan
x0 = [0 1.5; 0 -1.5; 1.5 0; -1 2.5; 2 2; -0.5 2; -0.2 1.1; 5 0; -1.5 2.5]';
% x0 = [0 1.5; 0 -1.5; 1.5 0; -1.5 0]';
TimeSpan = [0 1]; % The timespan should be small, for example [0 1].

bound.x1 = -3:0.1:3;
bound.x2 = bound.x1.^2 + 1;

%% ODE
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep',0.01);

for i = 1:size(x0,2)
    [T,X] = ode45(@fun_Matlab,TimeSpan,x0(:,i),options);
    t{i} = T;
    x{i} = X;
end

%% Plots
% figure
% plot(T,X(:,1),'r-')
% hold on
% plot(T,X(:,2),'b--')
% grid on
% legend('x_1','x_2')
% xlabel('time')
% ylabel('x_i')

figure
hold on
plot(bound.x1,bound.x2,'k--','linewidth',1.1)
for i = 1:size(x0,2)
    plot(x{i}(:,1),x{i}(:,2),'r-')
end
grid on
xlim([-2 2])
ylim([-1 3])
xlabel('$x_1$')
ylabel('$x_2$')