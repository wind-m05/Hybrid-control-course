clear all
close all

%% Calculation of Vectorfields in different modes
dx = 0.1;
% Right halve, x1 >= 0
[xp11,xp12] = meshgrid(-1.5:dx:1.5, -1.5:dx:1.5);
p11 = -2*xp11 + 4*xp12;
p12 = -2*xp12;
% Left halve, x1 < 0
[xp21,xp22] = meshgrid(-1.5:dx:1.5, -1.5:dx:1.5);
p21 = -4*xp21 - 7*xp22;
p22 = 4*xp21 - 5*xp22;

% Plotting of the vector fields
figure
quiver(xp11,xp12,p11,p12,'r-')
hold on    
quiver(xp21,xp22,p21,p22,'k-')
grid on
title('Vectorfield mode 2')
xlabel('x_1')
ylabel('x_2')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
xtest = linspace(-1.5,1.5);
plot(xtest,xtest.^2+1)
