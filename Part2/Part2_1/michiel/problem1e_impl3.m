close all; clc; clear 
A1 = [1.6 2; -4 -1.8];
A2 = [-1.8 4; -2 1.6];


%% grid search for theta
theta = 2:0.005:3;

%% Fix tau1 and tau2 
Phist = cell(length(theta),1);
for i = 1:length(theta)
yalmip clear
P = sdpvar(2,2,'symmetric');

tau1_l = 1.2;
tau1_h = 1.25;

tau2_l = 1;
tau2_h = theta(i);

% Implementation 1
% Abar1 = expm(A2*tau2_l) * expm(A1*tau1_l);
% Abar2 = expm(A2*tau2_l) * expm(A1*tau1_h);
% Abar3 = expm(A2*tau2_h) * expm(A1*tau1_l);
% Abar4 = expm(A2*tau2_h) * expm(A1*tau1_h);


% Implementation 2 -> lead to higher theta values, from 2.67 -> 2.78
Abar1 = expm(A2*tau2_h) * expm(A1*tau1_h) * expm(A2*tau2_l) * expm(A1*tau1_l);
Abar2 = expm(A2*tau2_l) * expm(A1*tau1_h) * expm(A2*tau2_h) * expm(A1*tau1_l);
Abar3 = expm(A2*tau2_l) * expm(A1*tau1_h) * expm(A2*tau2_l) * expm(A1*tau1_l);
Abar4 = expm(A2*tau2_h) * expm(A1*tau1_h) * expm(A2*tau2_h) * expm(A1*tau1_l);

% Implementation 3 Lead to higher theta values, from 2.78 -> 3, cant be
% correct, because you dont have this specific sequence periodically
Abar = Abar1*Abar2*Abar3*Abar4;

L1 = Abar'*P*Abar-P <=-1e9;
Lp = P >=1e-9;
L = L1+Lp;

opts = sdpsettings('solver','mosek','verbose',0);
optimize(L,[],opts);
% disp(diagnostics.info)
% if diagnostics.problem == 0
%     disp('Feasible')
% elseif diagnostics.problem == 1
%     disp('Infeasible')
% else
%     disp('Something else happened')
% end
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

%% Retrieve the feasible value range of theta
indx = find(defcheck == 1);
theta_val = theta(indx)