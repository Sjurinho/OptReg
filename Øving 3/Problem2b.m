clc;clear all;
x1 = linspace(0,4);
x2 = linspace(0,5);

[X, Y] = meshgrid(x1,x2);
f = (-3/2)*X - Y;



contour(X,Y,f);
xlabel("x1");
ylabel("x2");
title("Contour lines and constraints")
hold on;    
c1 = 15-3*x2;
c2 = 8-2*x1;
plot(c1,x2, "-k");
plot(x1,c2, "-k");
xlim([0,4]);
ylim([0,5]);

%%Find optimal solution
A = [2 1 1 0;
    1 3 0 1];
b = [8 15]';
c = [-3 -2 0 0]';
% A fesible starting point:
x0 = [0 0 8 15]'; % The starting point must also have four variables!
[x, fval, iterates] = simplex(c,A,b,x0,'report');

% Extract iterates in the space of x_1 and x_2:
iter_x1_x2 = iterates(1:2, :);

% Extract iterates as individual vectors (containing x_1 and x_2):
iter_1 = iter_x1_x2(:,1);
iter_2 = iter_x1_x2(:,2);
iter_3 = iter_x1_x2(:,3);