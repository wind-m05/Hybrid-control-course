clear
close all
clc
set(0,'defaulttextInterpreter','latex')
%% Parameters
Ap = [-0.2 0; 1 0];
Bp = [1; 0];
Cp = [1 0.1];
Dp = 1;



% Reference and disturbance are zero, so only A matrix in dynamics
r = 0;
d = 0;
a = (-3:0.1:3);
for i = 1:length(a)
yalmip clear
% Modes
A1 = [Ap Bp; -Cp a(i)];
A2 = [Ap -Bp; -Cp a(i)];

% S procedure
E_1_1 = [-0.1, -1, 0;
          0 , 0, 1];
E_1_2 = -E_1_1;
E_2_1 = [-0.1, -1, 0;
          0 , 0, -1];
E_2_2 = -E_2_1;

Z_1_12 = [-0.1,0,0;
      1  ,0,0;
      0  ,1,0];
Z_1_21 = Z_1_12;

Z_2_12 = [1 ,0 , 0;
      0 ,1 , 0;
      0 ,0 , 0]; 
Z_2_21 = Z_2_12;

% LMI variables
Pvar1 = sdpvar(3,3);
Pvar2 = sdpvar(3,3);

Uvar1 = sdpvar(2,2);
Uvar2 = sdpvar(2,2);

Wvar1 = sdpvar(2,2);
Wvar2 = sdpvar(2,2);

% LMIs
Lf1_1 = A1'*Pvar1 + Pvar1*A1 + E_1_1'*Uvar1*E_1_1 <= -1e-9;
Lf1_2 = A1'*Pvar1 + Pvar1*A1 + E_1_2'*Uvar1*E_1_2 <= -1e-9;
Lf2_1 = A2'*Pvar2 + Pvar2*A2 + E_2_1'*Uvar2*E_2_1 <= -1e-9;
Lf2_2 = A2'*Pvar2 + Pvar2*A2 + E_2_2'*Uvar2*E_2_2 <= -1e-9;

Lp1_1 = Pvar1 - E_1_1'*Wvar1*E_1_1 >= 1e-9;
Lp1_2 = Pvar1 - E_1_2'*Wvar1*E_1_2 >= 1e-9;
Lp2_1 = Pvar2 - E_2_1'*Wvar2*E_2_1 >= 1e-9;
Lp2_2 = Pvar2 - E_2_2'*Wvar2*E_2_2 >= 1e-9;

Lc_1_12 = Z_1_12'*(Pvar1 - Pvar2)*Z_1_12 == 0;
Lc_1_21 = Z_1_21'*(Pvar2 - Pvar1)*Z_1_21 == 0;
Lc_2_12 = Z_2_12'*(Pvar1 - Pvar2)*Z_2_12 == 0;
Lc_2_21 = Z_2_21'*(Pvar2 - Pvar1)*Z_2_21 == 0;

Lpos = [Uvar1(:)>=1e-9, Uvar2(:)>=1e-9, Wvar1(:)>=1e-9 , Wvar2(:)>=1e-9];
L = Lf1_1 + Lf1_2 + Lf2_1 + Lf2_2 + Lp1_1 + Lp1_2 + Lp2_1 + Lp2_2 + Lc_1_12 + Lc_1_21 + Lc_2_12 + Lc_2_21 + Lpos;

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

negdef_check = 0;
%% Should be negative definite

if all(eig(A1'*P1 + P1*A1 + E_1_1'*U1*E_1_1) < 0) && all(eig(A1'*P1 + P1*A1 + E_1_2'*U1*E_1_2) < 0) && all(eig(A2'*P2 + P2*A2 + E_2_1'*U2*E_2_1) <0) && all(eig(A2'*P2 + P2*A2 + E_2_2'*U2*E_2_2) < 0)
% eig(A1'*P1 + P1*A1 + E_1_1'*U1*E_1_1)
% eig(A1'*P1 + P1*A1 + E_1_2'*U1*E_1_2)
% eig(A2'*P2 + P2*A2 + E_2_1'*U2*E_2_1)
% eig(A2'*P2 + P2*A2 + E_2_2'*U2*E_2_2)
negdef_check = 1;
end

posdef_check = 0;
%% Should be positive definite
if all(eig(P1 - E_1_1'*W1*E_1_1) > 0) && all(eig(P1 - E_1_2'*W1*E_1_2) > 0) && all(eig(P2 - E_2_1'*W2*E_2_1) > 0) && all(eig(P2 - E_2_2'*W2*E_2_2) > 0)
% eig(P1 - E_1_1'*W1*E_1_1)
% eig(P1 - E_1_2'*W1*E_1_2)
% eig(P2 - E_2_1'*W2*E_2_1)
% eig(P2 - E_2_2'*W2*E_2_2)
posdef_check = 1;
end

%% Should be as close to zero as possible
Z_c1 = Z_1_12'*(P1 - P2)*Z_1_12;
Z_c2 = Z_1_21'*(P2 - P1)*Z_1_21;
Z_c3 = Z_2_12'*(P1 - P2)*Z_2_12;
Z_c4 = Z_2_21'*(P2 - P1)*Z_2_21;
Z_check{i} = [Z_c1 , Z_c2;
              Z_c3,  Z_c4];
all_check = 0;
if negdef_check && posdef_check
all_check = 1;
end
checker(i) = all_check;
end
indx = find(checker == 0);
a_optimal = a(indx(1));
