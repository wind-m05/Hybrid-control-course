clear
close all
clc
yalmip clear

% System dynamics:
A1 = [-7  4 6; 8 -47 -60; 0 36 45];
A2 = [ 7 16 8; 3   5  -1; 1  6 11];

% LMI variables:
Pvar = sdpvar(3,3); % symmetric by default

%% Search for a common quadratic Lyapunov function
% Lyapunov condition:
Lf1 = A1'*Pvar+Pvar*A1 <= -1e-1;
Lf2 = A2'*Pvar+Pvar*A2 <= -1e-1;

Lp = Pvar >= 1e-9;

% combine constraints into one object
% L = [Lf1,Lp]; %
L = [Lf1,Lf2,Lp]; 

% solve the LMI using SDPT3:
opts = sdpsettings('solver','mosek');
diagnostics = optimize(L,[],opts);
disp(diagnostics.info)
if diagnostics.problem == 0
 disp('Feasible')
elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

P = value(Pvar)
eigP = eig(P)
eigLf1 = eig(A1'*P+P*A1)
eigLf2 = eig(A2'*P+P*A2)

check(L)