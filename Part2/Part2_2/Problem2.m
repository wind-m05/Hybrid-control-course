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
Y12var = sdpvar(1,2);
Y21var = sdpvar(1,2);

% Conditions
Lp = Zvar >= 1e-9;
Lf12 = [Zvar (Zvar*A1' + Y12var'*B1'); (A1*Zvar + B1*Y12var) Zvar] >= 0;
Lf2 = [Zvar (Zvar*A2' + Y21var'*B2'); (A2*Zvar + B2*Y21var) Zvar] >= 0;
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

Y1 = value(Y12var);
Y2 = value(Y21var);
K1 = Y1*P;
K2 = Y2*P;
disp(eig(Z - (Z*A1' + Y1'*B1')*P*(A1*Z + B1*Y1)));
disp(eig(Z - (Z*A2' + Y2'*B2')*P*(A2*Z + B2*Y2)));

%% b Common QLF
pause
clear
close all
clc

A1 = [1/2 1/2; -1/4 -5/4];
A2 = [5/6 5/3; -2/3  1/6];
B1 = [ 0; 3];
B2 = [-2; 1];
C1 = [0 1];
C2 = [1/3 2/3];

nx = size(A1,1);

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
Y11var = sdpvar(2,1);
Y12var = sdpvar(2,1);
Y21var = sdpvar(2,1);
Y22var = sdpvar(2,1);

% Conditions
Lp1 = P1var >= 1e-9;
Lp2 = P2var >= 1e-9;
Lf11 = [P1var (A1'*P1var - C1'*Y11var'); (P1var*A1 - Y11var*C1) P1var] >= 0;
Lf12 = [P1var (A1'*P2var - C1'*Y12var'); (P2var*A1 - Y12var*C1) P2var] >= 0;
Lf21 = [P2var (A2'*P1var - C2'*Y21var'); (P1var*A2 - Y21var*C2) P1var] >= 0;
Lf22 = [P2var (A2'*P2var - C2'*Y22var'); (P2var*A2 - Y22var*C2) P2var] >= 0;
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

Y11 = value(Y11var);
Y12 = value(Y12var);
Y21 = value(Y21var);
Y22 = value(Y22var);
L1 = inv(P2)*Y12;
L2 = inv(P1)*Y21;
disp(eig(P1 - (A1' - C1'*L1')*inv(P1)*(A1 - L1*C1)));
disp(eig(P1 - (A1' - C1'*L1')*inv(P2)*(A1 - L1*C1)));
disp(eig(P2 - (A2' - C2'*L2')*inv(P1)*(A2 - L2*C2)));
disp(eig(P2 - (A2' - C2'*L2')*inv(P2)*(A2 - L2*C2)));
