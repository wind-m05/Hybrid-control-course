clear
close all
clc

%% Dynamics
A1 = [1/2 1/2; -1/4 -5/4];
A2 = [5/6 5/3; -2/3  1/6];
B1 = [ 0; 3];
B2 = [-2; 1];
C1 = [0 1];
C2 = [1/3 2/3];

nx = size(A1,1);

%% a
% LMI Variables
Zvar = sdpvar(2,2);
X12var = sdpvar(1,2);
X21var = sdpvar(1,2);

% Conditions
Lp = Zvar >= 1e-9;
Lf12 = [Zvar (Zvar*A1' + X12var'*B1'); (A1*Zvar + B1*X12var) Zvar] >= 0;
Lf2 = [Zvar (Zvar*A2' + X21var'*B2'); (A2*Zvar + B2*X21var) Zvar] >= 0;
L = [Lf12,Lf2,Lp]; 

% Solve
opts = sdpsettings('solver','sdpt3');
diagnostics = optimize(L,[],opts);
disp(diagnostics.info)
if diagnostics.problem == 0
 disp('Feasible')
elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

check(L)

pause(1)

% Results
Z = value(Zvar);
P = Z\eye(nx); % faster than inv(Z)
disp(eig(P))

Y1 = value(X12var);
Y2 = value(X21var);
K1 = Y1*P;
K2 = Y2*P;
disp(eig(Z - (Z*A1' + Y1'*B1')*P*(A1*Z + B1*Y1)));
disp(eig(Z - (Z*A2' + Y2'*B2')*P*(A2*Z + B2*Y2)));

%% b Common QLF
% observability
% W = sym('W',[2 2]);
% obsv  = A1'*W*A1 - W + C1'*C1 == 0;
% obsvCheck = solve(obsv,W);
% Wresult = [double(obsvCheck.W1_1) double(obsvCheck.W1_2); double(obsvCheck.W2_1) double(obsvCheck.W2_2)];
% eig(Wresult)
O1 = rref([C1; C1*A1]);
O2 = rref([C2; C2*A2]);

% LMI Variables
P1var = sdpvar(2,2);
P2var = sdpvar(2,2);
X11var = sdpvar(2,1);
X12var = sdpvar(2,1);
X21var = sdpvar(2,1);
X22var = sdpvar(2,1);

% Conditions
Lp1 = P1var >= 1e-9;
Lp2 = P2var >= 1e-9;
Lf11 = [P1var (A1'*P1var - C1'*X11var'); (P1var*A1 - X11var*C1) P1var] >= 0;
Lf12 = [P1var (A1'*P2var - C1'*X12var'); (P2var*A1 - X12var*C1) P2var] >= 0;
Lf21 = [P2var (A2'*P1var - C2'*X21var'); (P1var*A2 - X21var*C2) P1var] >= 0;
Lf22 = [P2var (A2'*P2var - C2'*X22var'); (P2var*A2 - X22var*C2) P2var] >= 0;
L = [Lf11,Lf12,Lf21,Lf22,Lp1,Lp2]; 

% Solve
opts = sdpsettings('solver','sdpt3');
diagnostics = optimize(L,[],opts);
disp(diagnostics.info)
if diagnostics.problem == 0
 disp('Feasible')
elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

check(L)

pause(1)

% Results
P1 = value(P1var);
disp(eig(P1))
P2 = value(P2var);
disp(eig(P2))

X11 = value(X11var);
X12 = value(X12var);
X21 = value(X21var);
X22 = value(X22var);
L1 = inv(P2)*X12;
L2 = inv(P1)*X21;
disp(eig(P1 - (A1' - C1'*L1')*inv(P1)*(A1 - L1*C1)));
disp(eig(P1 - (A1' - C1'*L1')*inv(P2)*(A1 - L1*C1)));
disp(eig(P2 - (A2' - C2'*L2')*inv(P1)*(A2 - L2*C2)));
disp(eig(P2 - (A2' - C2'*L2')*inv(P2)*(A2 - L2*C2)));

%% c
% Closed loop system with controller from (a) and observer from (b)
Acls1 = [(A1 + B1*K1) -B1*K1; zeros(nx,nx) (A1 - L1*C1)];
Acls2 = [(A2 + B2*K2) -B2*K2; zeros(nx,nx) (A2 - L2*C2)];