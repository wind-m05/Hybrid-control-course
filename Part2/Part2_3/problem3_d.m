clear
close all
clc
set(0,'defaulttextInterpreter','latex')
%% Parameters
Ap = [-0.2 0; 1 0];
Bp = [1; 0];
Cp = [1 0.1];
Dp = 1;

%% d
yalmip clear
% Reference and disturbance are zero, so only A matrix in dynamics
r = 0;
d = 0;

a = (-3:0.1:3);
gamma = (-1:0.1:1);
for i = 1:length(a)
%     for j = 1:length(gamma)
yalmip clear

% Modes
A1 = [Ap Bp; -Cp a(i)];
A2 = [Ap -Bp; -Cp a(i)];
% LMI variables
Pvar = sdpvar(3,3);
gamma = sdpvar(1,1);
% Auxillary matrices
Cbar = [Cp'*Cp, zeros(2,1);
        zeros(1,3)];
Bbar = [Bp;0];
% CL_aux = [A1'*Pvar + Pvar*A1 + Cbar , Pvar*Bbar;
%           Bbar'*Pvar                , -gamma(j)^2];
CL_aux1_var = [A1'*Pvar + Pvar*A1 + Cbar , Pvar*Bbar;
          Bbar'*Pvar                , -gamma^2];

CL_aux2_var = [A2'*Pvar + Pvar*A2 + Cbar , Pvar*Bbar;
          Bbar'*Pvar                , -gamma^2];

% LMIs
Lf1 = CL_aux1_var <= -1e-9;
Lf2 = CL_aux2_var <= -1e-9;
Lp = Pvar >= 1e-9;
L = Lf1 + Lf2 + Lp;

opts = sdpsettings('solver','mosek','verbose',1);
optimize(L,[],opts);
% To manually check constraints
P = value(Pvar);

CL_aux1 = [A1'*P + P*A1 + Cbar , P*Bbar;
          Bbar'*P               , -gamma(j)^2];

CL_aux2 = [A2'*P + P*A2 + Cbar , P*Bbar;
          Bbar'*P                , -gamma(j)^2];
posdef_check = 0;
if all(eig(P) > 0)
posdef_check = 1;
end
negdef_check = 0;
if all(eig(CL_aux1) < 0 ) && all(eig(CL_aux2) < 0 )
negdef_check  = 1;
end
%% Should be negative definite

% if all(eig(A1'*P1 + P1*A1 + E_1_1'*U11*E_1_1) < 0) && all(eig(A1'*P1 + P1*A1 + E_1_2'*U12*E_1_2) < 0) && all(eig(A2'*P2 + P2*A2 + E_2_1'*U21*E_2_1) <0) && all(eig(A2'*P2 + P2*A2 + E_2_2'*U22*E_2_2) < 0)
% eig(A1'*P1 + P1*A1 + E_1_1'*U1*E_1_1)
% eig(A1'*P1 + P1*A1 + E_1_2'*U1*E_1_2)
% eig(A2'*P2 + P2*A2 + E_2_1'*U2*E_2_1)
% eig(A2'*P2 + P2*A2 + E_2_2'*U2*E_2_2)
% negdef_check = 1;
% end

% all_check = 0;
% if negdef_check && posdef_check
% all_check = 1;
% end
% checker(i,j) = all_check;
%     end
end
% indx = find(checker == 0);
% a_optimal = a(indx(1));