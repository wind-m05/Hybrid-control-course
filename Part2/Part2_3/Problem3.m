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

%% b
r = 1;
d = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,x] = ode45(@(t,y) fun3(t,y,Ap,Bp,Cp,Dp,a,r,d), tspan, x0, options);
[t_i,x_i] = ode45(@(t,y) int3(t,y,Ap,Bp,Cp,Dp,a,r,d), tspan, x0, options);

figure
subplot(3,1,1)
plot(t,x(:,1)); grid on; hold on
plot(t_i,x_i(:,1))
ylabel('$x_{p,1}$')
subplot(3,1,2)
plot(t,x(:,2)); grid on; hold on
plot(t_i,x_i(:,2))
ylabel('$x_{p,2}$')
subplot(3,1,3)
plot(t,x(:,3)); grid on; hold on
plot(t_i,x_i(:,3))
ylabel('$x_{1}$')
xlabel('Time $t$ [s]')