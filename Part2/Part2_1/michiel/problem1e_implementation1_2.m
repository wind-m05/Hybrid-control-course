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


% Theta = 2.720 if you only include Abar1 -> Abar4, 
% theta = 2.595 if you include Abar5 -> Abar8
Abar1 = expm(A2*tau2_h) * expm(A1*tau1_l) * expm(A2*tau2_l) * expm(A1*tau1_h);
Abar2 = expm(A2*tau2_l) * expm(A1*tau1_l) * expm(A2*tau2_h) * expm(A1*tau1_h);
Abar3 = expm(A2*tau2_l) * expm(A1*tau1_l) * expm(A2*tau2_l) * expm(A1*tau1_h);
Abar4 = expm(A2*tau2_h) * expm(A1*tau1_l) * expm(A2*tau2_h) * expm(A1*tau1_h);

Abar5 = expm(A1*tau1_l) * expm(A2*tau2_l) * expm(A1*tau1_h) * expm(A2*tau2_h);
Abar6 = expm(A1*tau1_l) * expm(A2*tau2_h) * expm(A1*tau1_h) * expm(A2*tau2_l); 
Abar7 = expm(A1*tau1_l) * expm(A2*tau2_l) * expm(A1*tau1_h) * expm(A2*tau2_l); 
Abar8 = expm(A1*tau1_l) * expm(A2*tau2_h) * expm(A1*tau1_h) * expm(A2*tau2_h);


L1 = Abar1'*P*Abar1-P <=-1e9;
L2 = Abar2'*P*Abar2-P <=-1e9;
L3 = Abar3'*P*Abar3-P <=-1e9;
L4 = Abar4'*P*Abar4-P <=-1e9;

L5 = Abar5'*P*Abar5-P <=-1e9;
L6 = Abar6'*P*Abar6-P <=-1e9;
L7 = Abar7'*P*Abar7-P <=-1e9;
L8 = Abar8'*P*Abar8-P <=-1e9;


Lp = P >=1e-9;
L = L1+L2+L3+L4+L5+L6+L7+L8+Lp;

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