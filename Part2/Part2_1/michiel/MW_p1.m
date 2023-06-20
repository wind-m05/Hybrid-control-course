clear
close all
clc

%% Dynamics
A1 = [1.6 2; -4 -1.8];
A2 = [-1.8 4; -2 1.6];
nx = size(A1,1);

% a
alpha = 0:0.01:1;
As = zeros(nx,nx,length(alpha));
for i = 1:length(alpha)
    As(:,:,i) = alpha(i)*A1 + (1-alpha(i))*A2;
end

figure
hold on
for i = 1:length(alpha)
    plot(eig(As(:,:,i)),'ro')
end
xlim([-1,1])
ylim([-3,3])
grid on

%% c
A1 = [1.6 2; -4 -1.8];
A2 = [-1.8 4; -2 1.6];
nx = size(A1,1);
tau1 = 0:0.001:1.5;
tau2 = tau1;

% Check Hurwitz 
delta_check = zeros(length(tau1),length(tau2),1);
delta_mat = zeros(length(tau1),length(tau2),1);
for p = 1:length(tau1)
    for q = 1:length(tau2)
        delta =  tau1(p)/(tau1(p)+tau2(q));
        delta_mat(p,q) = delta;
        A_cc = delta*A1+(1-delta)*A2;
        if eig(A_cc(1)) < 0 && eig(A_cc(2)) < 0 
        delta_check(p,q) = 1;
        else
        delta_check(p,q) = 0;
        end
    end
end

schur_check = zeros(length(tau1),length(tau2),1);
% Check Schur
for p = 1:length(tau1)
    for q = 1:length(tau2)

        A_sh = exp(A2*(1-delta_mat(p,q))*(tau1(p)+tau2(q))) * exp(A1*(delta_mat(p,q)*(tau1(p)+tau2(q))));
        if abs(eig(A_sh(1))) < 1 && abs((A_sh(2))) < 1 
        schur_check(p,q) = 1;
        else
        schur_check(p,q) = 0;
        end
    end
end
final_check = schur_check.*delta_check;

%% d