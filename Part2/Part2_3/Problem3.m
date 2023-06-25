clear
close all
clc
set(0,'defaulttextInterpreter','latex')

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

% Plot
figure
subplot(4,1,1)
plot(t,x(:,1)); grid on; hold on
plot(t_i,x_i(:,1))
xline(e_trigger,'--')
xlim([0 25])
% ylim([-1.5 2.5])
ylabel('$x_{p,1}(t)$')
legend('SPAN','LTI integrator','interpreter','latex')

subplot(4,1,2)
plot(t,x(:,2)); grid on; hold on
plot(t_i,x_i(:,2))
xline(e_trigger,'--')
xlim([0 25])
% ylim([6 11])
ylabel('$x_{p,2}(t)$')

subplot(4,1,3)
plot(t,x(:,3)); grid on; hold on
plot(t_i,x_i(:,3))
xline(e_trigger,'--')
xlim([0 25])
% ylim([-1.5 1.5])
ylabel('$x_{1}(t)$')

subplot(4,1,4)
plot(t,e); grid on; hold on
plot(t_i,e_i)
xline(e_trigger,'--')
xlim([0 25])
% ylim([-2 2])
ylabel('$e(t)$')
xlabel('Time $t$ [s]')

%% c
yalmip clear
% Reference and disturbance are zero, so only A matrix in dynamics
r = 0;
d = 0;

% [x1, xp2] = meshgrid(-5:0.1:5);
% xp1 = -0.1*xp2;
% 
% figure
% hold on
% surf(xp1, xp2, x1,'EdgeColor','none','Facecolor',[0.9 0.1 0.1],'FaceAlpha',0.5)
% % scatter3(x(:,1),x(:,2),x(:,3))
% grid on
% xlabel('xp1')
% ylabel('xp2')
% zlabel('x1')
% xlim([-3 3])
% ylim([-3 3])
% zlim([-3 3])

% Varying a NEEDS TO BE ADDED

% Modes
A1 = [Ap Bp; -Cp a];
A2 = [Ap -Bp; -Cp a];

% S procedure
E2 = [-1 -0.1 0; 0 0 1];
E1 = -E2;
Z12 = [1; -10; 0];
Z21 = Z12;

% LMI variables
Pvar1 = sdpvar(3,3);
Pvar2 = sdpvar(3,3);
Uvar1 = sdpvar(2,2);
Uvar2 = sdpvar(2,2);
Wvar1 = sdpvar(2,2);
Wvar2 = sdpvar(2,2);

% LMIs
Lf1 = A1'*Pvar1 + Pvar1*A1 + E1'*Uvar1*E1 <= -1e-9;
Lf2 = A2'*Pvar2 + Pvar2*A2 + E2'*Uvar2*E2 <= -1e-9;
Lp1 = Pvar1 - E1'*Wvar1*E1 >= 1e-9;
Lp2 = Pvar2 - E2'*Wvar2*E2 >= 1e-9;
Lc12 = Z12'*(Pvar1 - Pvar2)*Z12 == 0;
% Lc12g = Z12'*(Pvar1 - Pvar2)*Z12 >= -1e-7;
% Lc12l = Z12'*(Pvar1 - Pvar2)*Z12 <= 1e-7;
Lc21 = Z21'*(Pvar2 - Pvar1)*Z21 == 0;
% Lc21g = Z21'*(Pvar2 - Pvar1)*Z21 >= -1e-7;
% Lc21l = Z21'*(Pvar2 - Pvar1)*Z21 <= 1e-7;
Lpos = [Uvar1(:)>=1e-9, Uvar2(:)>=1e-9, Wvar1(:)>=1e-9, Wvar2(:)>=1e-9];
L = Lf1 + Lf2 + Lp1 + Lp2 + Lc12 + Lc21 + Lpos;
% L = Lf1 + Lf2 + Lp1 + Lp2 + Lc12g + Lc12l + Lc21g + Lc21l + Lpos;

opts = sdpsettings('solver','mosek','verbose',2);
optimize(L,[],opts);
% disp(diagnostics.info)
% if diagnostics.problem == 0
%  disp('Feasible')
% elseif diagnostics.problem == 1
%  disp('Infeasible')
% else
%  disp('Something else happened')
% end

% To manually check constraints
P1 = value(Pvar1);
P2 = value(Pvar2);
U1 = value(Uvar1);
U2 = value(Uvar2);
W1 = value(Wvar1);
W2 = value(Wvar2);

eig(A1'*P1 + P1*A1 + E1'*U1*E1)
eig(A2'*P2 + P2*A2 + E2'*U2*E2)
eig(P1 - E1'*W1*E1)
eig(P2 - E2'*W2*E2)
Z12'*(P1 - P2)*Z12
Z21'*(P2 - P1)*Z21