clear
close all
clc

%% Search for a Piecewise quadratic one with S-procedure: P1, P2, U1, U2, W1, W2
A1 = [ 2  4; 0 -1];
A2 = [-4 -7; 4 -5];

% regions
E1=[-2 1; 2 1];
E3=-E1;
E2=[2 -1; 2 1];
E4=-E2;

% continuity
Z12 = [ 1; 2];
Z14 = [-1; 2];
Z21 = Z12;
Z34 = Z12;
Z43 = Z12;
Z41 = Z14;
Z23 = Z14;
Z32 = Z14;

% LMI variables:
Pvar1 = sdpvar(2,2);
Pvar2 = sdpvar(2,2);
Pvar3 = sdpvar(2,2);
Pvar4 = sdpvar(2,2);

Uvar1=sdpvar(2,2);
Uvar2=sdpvar(2,2);
Uvar3=sdpvar(2,2);
Uvar4=sdpvar(2,2);

Wvar1=sdpvar(2,2);
Wvar2=sdpvar(2,2);
Wvar3=sdpvar(2,2);
Wvar4=sdpvar(2,2);

% Lyapunov condition:
Lf1 = A1'*Pvar1 + Pvar1*A1 + E1'*Uvar1*E1 <= -1e-9;
Lf2 = A2'*Pvar2 + Pvar2*A2 + E2'*Uvar2*E2 <= -1e-9;
Lf3 = A1'*Pvar3 + Pvar3*A1 + E3'*Uvar3*E3 <= -1e-9;
Lf4 = A2'*Pvar4 + Pvar4*A2 + E4'*Uvar4*E4 <= -1e-9;

Lp1 = Pvar1 - E1'*Wvar1*E1 >= 1e-9;
Lp2 = Pvar2 - E2'*Wvar2*E2 >= 1e-9;
Lp3 = Pvar3 - E3'*Wvar3*E3 >= 1e-9;
Lp4 = Pvar4 - E4'*Wvar4*E4 >= 1e-9;

Lc12 = Z12'*(Pvar1 - Pvar2)*Z12 == 0;
Lc14 = Z14'*(Pvar1 - Pvar4)*Z14 == 0;
Lc21 = Z21'*(Pvar2 - Pvar1)*Z21 == 0;
Lc23 = Z23'*(Pvar2 - Pvar3)*Z23 == 0;
Lc32 = Z32'*(Pvar3 - Pvar2)*Z32 == 0;
Lc34 = Z34'*(Pvar3 - Pvar4)*Z34 == 0;
Lc41 = Z41'*(Pvar4 - Pvar1)*Z41 == 0;
Lc43 = Z43'*(Pvar4 - Pvar3)*Z43 == 0;

Lpos = [Uvar1(:)>=0, Uvar2(:)>=0, Uvar3(:)>=0, Uvar4(:)>=0, Wvar1(:)>=0, Wvar2(:)>=0, Wvar3(:)>=0, Wvar4(:)>=0];

L = Lf1 + Lf2 + Lf3 + Lf4 + ...
    Lp1 + Lp2 + Lp3 + Lp4 +...
    Lc12 + Lc14 + Lc21 + Lc23 + Lc32 + Lc34 + Lc41 + Lc43 +...
    Lpos;

%solve the LMI using SDPT3:
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