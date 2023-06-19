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

% figure
% hold on
% for i = 1:length(alpha)
%     plot(eig(As(:,:,i)),'ro')
% end
% xlim([-1,1])
% ylim([-3,3])
% grid on

% c
tau1 = 0:0.1:1.5;
tau2 = tau1;

stabAk = zeros(length(tau1),length(tau2),1);
for p = 1:length(tau1)
    for q = 1:length(tau2)
        Ak = exp(A2*tau2(q)) * exp(A1*tau1(p))
        Ak = (eye(nx) + A2*tau2(q))*(eye(nx) + A1*tau1(p))
        eigAk = eig(Ak);
        if abs(eigAk(1)) < 1 && abs(eigAk(2)) < 1
            stabAk(p,q,:) = 1; % Schur stable
        else
            stabAk(p,q,:) = 0; % unstable
        end
    end
end