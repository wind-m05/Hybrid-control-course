clear
close all
clc
set(0,'defaulttextInterpreter','latex')

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
tau1 = 0:0.01:15;
tau2 = 0:0.01:15;

stabAk = zeros(length(tau1),length(tau2),1);
for p = 1:length(tau1)
    for q = 1:length(tau2)
        Ak = expm(A2*tau2(q)) * expm(A1*tau1(p));
%         Ak = (eye(nx) + A2*tau2(q))*(eye(nx) + A1*tau1(p))
        eigAk = eig(Ak);
        if abs(eigAk(1)) < 1 && abs(eigAk(2)) < 1
            stabAk(p,q,:) = 1; % Schur stable
        else
            stabAk(p,q,:) = 0; 
            % unstable
        end
    end
end

%% Plot
map = [0.9 0.1 0.1; 0.05 0.7 0.05];

[Tau1, Tau2] = meshgrid(tau1, tau2);
figure
scatter(Tau1(:),Tau2(:),10,stabAk(:),'filled')
colormap(map)
colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',{'','Unstable','','Stable'.''},'TickLabelInterpreter','latex')
xlabel('Time in mode 1 $\tau_1$ [time units]')
ylabel('Time in mode 2 $\tau_2$ [time units]')
% xlim([0 1.5])
% ylim([0 1.5])

