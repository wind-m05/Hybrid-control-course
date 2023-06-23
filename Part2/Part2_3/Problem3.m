clear
close all
clc

%% Parameters
Ap = [-0.2 0; 1 0];
Bp = [1; 0];
Cp = [1 0.1];
Dp = 1;
a = 0;

tspan = [0 50];
x0 = [0 0 0]';
% x0 = [2 5 1]';

%% b
r = 1;
d = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,x] = ode45(@(t,y) fun3(t,y,Ap,Bp,Cp,Dp,a,r,d), tspan, x0, options);
[t_i,x_i] = ode45(@(t,y) int3(t,y,Ap,Bp,Cp,Dp,a,r,d), tspan, x0, options);

e = zeros(length(x),1);
e_i = zeros(length(x_i),1);
for i = 1:length(x)
    e(i) = -Cp*[x(i,1); x(i,2)] + Dp*r;
end
for i = 1:length(x_i)
    e_i(i) = -Cp*[x_i(i,1); x_i(i,2)] + Dp*r;
end

e_trigger = t(e<1e-3);
e_trigger = e_trigger(1);

%%
figure
subplot(4,1,1)
plot(t,x(:,1)); grid on; hold on
plot(t_i,x_i(:,1))
xline(e_trigger,'--')
xlim([0 25])
ylim([-1.5 2.5])
ylabel('$x_{p,1}(t)$')
legend('SPAN','LTI integrator','interpreter','latex')

subplot(4,1,2)
plot(t,x(:,2)); grid on; hold on
plot(t_i,x_i(:,2))
xline(e_trigger,'--')
xlim([0 25])
ylabel('$x_{p,2}(t)$')

subplot(4,1,3)
plot(t,x(:,3)); grid on; hold on
plot(t_i,x_i(:,3))
xline(e_trigger,'--')
xlim([0 25])
ylim([-1.5 1.5])
ylabel('$x_{1}(t)$')

subplot(4,1,4)
plot(t,e); grid on; hold on
plot(t_i,e_i)
xline(e_trigger,'--')
xlim([0 25])
ylabel('$e(t)$')
xlabel('Time $t$ [s]')
