close all;clc;clear

A1 = [-7 ,4   ,6;
       8 ,-47 ,-60;
       0 ,36  ,45];
A2 = [-4/3,13/3,6;
       -1/3,-8/3,-1;
        0,0,-2];

%% Check how the decay works
t = 0:0.01:10;
c = [1 5];
x0 = 10;
lambda = [1 5];
figure
for i = 1:length(lambda)
y = c(i)*exp(-lambda(i)*t)*norm(x0);
plot(t,y)
hold on
end
