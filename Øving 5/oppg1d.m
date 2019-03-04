%Getting constants
clear all;
run Constants.m

z = zeros(N*nx + N*mu,1);

%Minimization problem
KKT_mat = [G -Aeq';
    Aeq zeros(N*nx)];
KKT_vec = [zeros(1, N*nx + N*mu) beq]';

lambda = zeros(size(KKT_mat,1) - size(z,1),1);

z_and_lambda = KKT_mat^(-1)*KKT_vec;

z = z_and_lambda(1:N*nx+N*mu);

lambda = z_and_lambda(N*nx+N*mu+1:210);

%Extract solution
y = [x0(3) z(nx:3:N*nx)']';
u = z(N*nx+1:N*nx+N*mu);

%Plots
plot(0:N, y, '-o')
hold on;
plot(0:N-1, u, '-o')
title('Open loop optimized system')
xlabel('timestep[n]')
legend({'$y$', '$u$'}, 'Interpreter', 'Latex', 'FontSize', 14);
hold off;

