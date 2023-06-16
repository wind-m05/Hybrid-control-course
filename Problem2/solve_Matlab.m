clear
close all
clc

set(0,'defaulttextInterpreter','latex')

%% Dynamics
A1 = [2 4; 0 -1];
A2 = [-4 -7; 4 -5];

%% Initial conditions and timespan
% x0 = [0 1.5; 0 -1.5; 1.5 0; -1 2.5; 2 2; -0.5 2; -0.2 1.1; 5 0; -1.5 2.5]';
x0 = [1 2]';
TimeSpan = [0 1]; % The timespan should be small, for example [0 1].

bound.x1 = -3:0.1:3;
bound.x2 = bound.x1.^2 + 1;

%% ODE
options = odeset('RelTol',1e-1,'AbsTol',1e-1,'MaxStep',0.01);

for i = 1:size(x0,2)
    [T,X] = ode45(@fun_Matlab_lambda,TimeSpan,x0(:,i),options);
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
legend('Switching plane','Solutions','interpreter','latex')

%% Plot with regions
x1_plot = -2:0.01:2;

figure
hold on
plot(bound.x1,bound.x2,'k--','linewidth',1.1)
plot(x1_plot,2*x1_plot,'r')
plot(x1_plot,-2*x1_plot,'r')
text(-0.15,2.25,'1','FontSize',14,'EdgeColor','k')
text( 1.25,0.25,'2','FontSize',14,'EdgeColor','k')
text(-0.15,-0.75,'3','FontSize',14,'EdgeColor','k')
text(-1.25,0.25,'4','FontSize',14,'EdgeColor','k')
grid on
xlim([-2 2])
ylim([-1 3])
xlabel('$x_1$')
ylabel('$x_2$')
legend('Switching plane','Region boundaries','interpreter','latex')