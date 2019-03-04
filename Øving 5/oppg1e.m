%Getting constants
clear all;
run Constants.m

%[z, fval, exitflag, opts] = quadprog(G,[], [], [], Aeq, beq, [], [], x0); %1e
[z, fval, exitflag, opts] = quadprog(G,[], [], [], Aeq, beq, lb, ub, x0); %1f
%Extract solution

y = [x0(3) z(nx:3:N*nx)']';
u = z(N*nx+1:N*nx+N*mu);

disp(opts);

%Plots
figure(2)
plot(0:N, y, '-o')
hold on;
plot(0:N-1, u, '-o')
title('Open loop optimized system')
xlabel('timestep[n]')
legend({'$y$', '$u$'}, 'Interpreter', 'Latex', 'FontSize', 14);
%hold off;
