close all
clear all
A = [1 2;
     3 5];
R = [0.6 0.3;
     1   0.2];

Pvar = sdpvar(2,2);

Q = 0.1*eye(2); % Matrix related to the decrease rate
D = 0.1*eye(2);
% Define the constraints
Lp = Pvar >= 1e-9;
% Lu = Pvar <= 100;    % positive definiteness of  P
Lf = A'*Pvar+Pvar*A <= 0;
Lj = R'*Pvar*R-Pvar <= 0;
L = [Lp,Lf,Lj];                    % combine LMIs
% We use Pvar>= 1e-9 instead of Pvar>0, as strict inequalities do not make
% much sense in continuous numerical optimization.

% Solving the LMI conditions with SDPT3: optimize() creates the correct
% arguments and calls the SDPT3 solver to obtain the solution.
opts = sdpsettings('solver','sedumi ');
diagnostics = optimize(L,[],opts);
disp(diagnostics.info)
if diagnostics.problem == 0
    disp('Feasible')
elseif diagnostics.problem == 1
    disp('Infeasible')
else
    disp('Something else happened')
end
pause;

% value() converts the solution Pvar into a numeric matrix
P = value(Pvar)
eigP = eig(P)
disp('Eigenvalues are positive so P is ok, so the system is stable.')