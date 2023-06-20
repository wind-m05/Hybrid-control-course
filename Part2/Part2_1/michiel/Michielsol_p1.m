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
plot(eig(As(:,:,1)),'bo')
hold on
plot(eig(As(:,:,end)),'go')
 plot(eig(As(:,:,2)),'ro')
L = legend('A1','A2','$\alpha*A1 + (1-\alpha)*A2$','Interpreter','latex');
L.AutoUpdate = 'off'; 
hold on
for i = 3:length(alpha)-1
    plot(eig(As(:,:,i)),'ro')
    hold on
end
xline(0)
yline(0)
xlim([-1,1])
ylim([-3,3])
xlabel('$\sigma$','Interpreter','latex')
ylabel('$j\omega$','Interpreter','latex')
grid on



%% b
clear all
A1 = [-1,10;
      -100,-1];
A2 = [-1,100;
      -10,-1];
nx = size(A1,1);
alpha = 0:0.01:1;
As = zeros(nx,nx,length(alpha));
for i = 1:length(alpha)
    As(:,:,i) = alpha(i)*A1 + (1-alpha(i))*A2;
end

figure
plot(eig(As(:,:,1)),'bo')
hold on
plot(eig(As(:,:,end)),'go')
 plot(eig(As(:,:,2)),'ro')
L = legend('A1','A2','$\alpha*A1 + (1-\alpha)*A2$','Interpreter','latex');
L.AutoUpdate = 'off'; 
hold on
for i = 3:length(alpha)-1
    plot(eig(As(:,:,i)),'ro')
    hold on
end
xline(0)
yline(0)
xlim([-2,1])
ylim([-60,60])
xlabel('$\sigma$','Interpreter','latex')
ylabel('$j\omega$','Interpreter','latex')
grid on



%% c
clear all
A1 = [1.6 2; -4 -1.8];
A2 = [-1.8 4; -2 1.6];
nx = size(A1,1);
tau1 = 0:0.01:1.5;
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


% Check Schur
% epsilon = 0:0.1:(max(tau1)+max(tau2));
schur_check = zeros(length(tau1),length(tau2));
A_sh_eig = zeros(length(tau1),length(tau2));
for p = 1:length(tau1)
    for q = 1:length(tau2)
        A_sh = expm(A2*(1-delta_mat(p,q))*(tau1(p)+tau2(q))) * expm(A1*(delta_mat(p,q)*(tau1(p)+tau2(q))));
%         A_sh_eig(p,q,i) = eig(exp(A2*(1-delta_mat(p,q))*epsilon(i)) * exp(A1*(delta_mat(p,q)*epsilon(i))));
            if abs(eig(A_sh(1))) < 1 && abs((A_sh(2))) < 1 
            schur_check(p,q) = 1;
            else
            schur_check(p,q) = 0;
            end
    end
end
final_check = schur_check.*delta_check;
figure
imagesc(final_check)
colorbar
%% Test principle from the lecture notes similar to c 

A1 = [-0.5,1;
      100,-1];
A2 = [-1,-100;
      -0.5,-1];

% Take convex combination for switching time 1/2 epsilon

delta = 0.5;
A_cc = delta*A1 + (1-delta) *A2;
epsilon = 0:0.001:0.05;
A_sh = cell(length(epsilon),1);
for i = 1:length(epsilon)
    A_sh{i} = abs(eig(expm(1/2*epsilon(i)*A2)*expm(1/2*epsilon(i)*A1)));
plot(epsilon(i),A_sh{i})
hold on
end
