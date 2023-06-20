close all; clc; clear 
A1 = [1.6 2; -4 -1.8];
A2 = [-1.8 4; -2 1.6];
P = sdpvar(2,2);
Lp = P >= 0;
%% grid search for theta
theta = 1.5:0.01:2;

%% Fix tau1 and tau2 
Phist = cell(length(theta),1);
for i = 1:length(theta)
tau1_l = 0.15;
tau1_h = 0.2;
tau2_l = 1.5;
tau2_h = theta(i);
Abarll = expm(A2*tau2_l) * expm(A1*tau1_l);
Abarlh = expm(A2*tau2_l) * expm(A1*tau1_h);
Abarhl = expm(A2*tau2_h) * expm(A1*tau1_l);
Abarhh = expm(A2*tau2_h) * expm(A1*tau1_h);

Lll = Abarll'*P*Abarll-P <=0;
Llh = Abarlh'*P*Abarlh-P <=0;
Lhl = Abarlh'*P*Abarhl-P <=0;
Lhh = Abarhh'*P*Abarhh-P <=0;
Lp = P >=1e-9;
L = Lll+Llh+Lhl+Lhh+Lp;

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
Phist{i} = value(P);
end

defcheck = zeros(length(theta),1);
for i = 1:length(theta)
    P = Phist{i};
    eigP = eig(P);
    if eigP(1) > 1 && eigP(2) > 0
    defcheck(i) = 1;
    else
    defcheck(i) = 0;
    end
end