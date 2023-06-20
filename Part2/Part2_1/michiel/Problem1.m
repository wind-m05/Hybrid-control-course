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
A2 = [-1.8 4; -2 3];
nx = size(A1,1);


tau1 = 0:0.1:1.5;
tau2 = tau1;

stabAk = zeros(length(tau1),length(tau2),1);
for p = 1:length(tau1)
    for q = 1:length(tau2)
        Ak = expm(A2*tau2(q))* expm(A1*tau2(p));
%         Ak = (eye(nx) + A2*tau2(q))*(eye(nx) + A1*tau1(p));
        eigAk = eig(Ak);
        
        if abs(eigAk(1)) < 1 && abs(eigAk(2)) < 1
            stabAk(p,q,:) = 1; % Schur stable
        else
            stabAk(p,q,:) = 0; % unstable
        end
    end
end
figure
imagesc(stabAk)
colorbar


stabAk_rev = zeros(length(tau1),length(tau2),1);
for p = 1:length(tau1)
    for q = 1:length(tau2)
        Ak =  expm(A1*tau2(p)) * expm(A2*tau2(q));
%         Ak = (eye(nx) + A2*tau2(q))*(eye(nx) + A1*tau1(p));
        eigAk = eig(Ak);
        
        if abs(eigAk(1)) < 1 && abs(eigAk(2)) < 1
            stabAk_rev(p,q,:) = 1; % Schur stable
        else
            stabAk_rev(p,q,:) = 0; % unstable
        end
    end
end
figure
imagesc(stabAk_rev)
colorbar
isequal(stabAk,stabAk_rev)


%% Other implementation for c

tau1 = 0:0.005:1.5;
tau2 = tau1;

stabAk = zeros(length(tau1),length(tau2),1);
for p = 1:length(tau1)
    for q = 1:length(tau2)
%         delta = tau1(p)/(tau1(p)+tau2(q));
%         epsilon = tau1(p)+tau2(q);
        Ak = expm((tau1(p)*A1+tau2(q)*A2));
%         Ak = (eye(nx) + A2*tau2(q))*(eye(nx) + A1*tau1(p));
        eigAk = eig(Ak);
        
        if abs(eigAk(1)) < 1 && abs(eigAk(2)) < 1
            stabAk(p,q,:) = 1; % Schur stable
        else
            stabAk(p,q,:) = 0; % unstable
        end
        eigAk_hist{p,q} = eigAk;
    end
end
figure
imagesc(stabAk)
colorbar

