clear
close all
clc
set(0,'defaulttextInterpreter','latex')
%% Parameters
Ap = [-0.2 0; 1 0];
Bp = [1; 0];
Cp = [1 0.1];
Dp = 1;


yalmip clear
% Reference and disturbance are zero, so only A matrix in dynamics
r = 0;
d = 0;
a = (-3:0.1:3);
for i = 1:length(a)
% Modes
A1 = [Ap Bp; -Cp a(i)];
A2 = [Ap -Bp; -Cp a(i)];

% S procedure
% E_1_1 = [-0.1, -1, 0;
%           0 , 0, 1];
% E_1_2 = -E_1_1;
% E_2_1 = [-0.1, -1, 0;
%           0 , 0, -1];
% E_2_2 = -E_2_1;
% 
% Z_1_12 = [0.1,0,0;
%       1  ,0,0;
%       0  ,1,0];
% Z_1_21 = Z_1_12;
% 
% Z_2_12 = [1 ,0 , 0;
%       0 ,1 , 0;
%       0 ,0 , 0]; 
% Z_2_21 = Z_2_12;

% LMI variables
Pvar1 = sdpvar(3,3);
Pvar2 = sdpvar(3,3);

% Uvar1 = sdpvar(2,2);
% Uvar2 = sdpvar(2,2);
% 
% Wvar1 = sdpvar(2,2);
% Wvar2 = sdpvar(2,2);

% LMIs
Lf1_1 = A1'*Pvar1 + Pvar1*A1  <= -1e-9;
Lf2_1 = A2'*Pvar2 + Pvar2*A2  <= -1e-9;
LP1 = Pvar1 >=1e-9;
LP2 = Pvar2 >=1e-9;
L = Lf1_1 + Lf2_1 + LP1 + LP2;

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

negdef_check = 0;
%% Should be negative definite

if all(eig(A1'*P1 + P1*A1) < 0) && all(eig(A2'*P2 + P2*A2) < 0)
% eig(A1'*P1 + P1*A1 + E_1_1'*U1*E_1_1)
% eig(A1'*P1 + P1*A1 + E_1_2'*U1*E_1_2)
% eig(A2'*P2 + P2*A2 + E_2_1'*U2*E_2_1)
% eig(A2'*P2 + P2*A2 + E_2_2'*U2*E_2_2)
negdef_check = 1;
end

posdef_check = 0;
%% Should be positive definite
if all(eig(P1) > 0) && all(eig(P2) > 0)
posdef_check = 1;
end

%% Should be as close to zero as possible
all_check = 0;
if negdef_check && posdef_check
all_check = 1;
end
checker(i) = all_check;
end