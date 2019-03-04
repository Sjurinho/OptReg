%Get constants
clear all;
run Constants.m

%Full state information => C = [1 1 1] => y_t = x_t
%function x_(t+1) = h(x) = Ax_t + Bu_t
h = @(x,u) A*x+B*u;
%Start condition
Y = [x0];
U = [0];
y = x0;
for n = 0:N-1
    %%Solve QP-problem
    beq(1:nx) = A*y;
    z = quadprog(G,[], [], [], Aeq, beq, lb, ub);
    %Extract solution
    %y = z(1:nx); %NB, only one step saved (x1, x2, x3)
    u = z(N*nx + 1); %Same goes for input
    
    %%Set preconditions for next iteration
    y = h(y,u); %Next starting condition for quadprog (x_(t+1))

    
    Y = [Y y]; %Save answers for each iteration
    U = [U u]; %Same for inputs
    
end

%Extract y3 and u
y_mpc = 1:N+1;
for i = 1:N+1
    y_mpc(i) = C*Y(:,i);
end
u_mpc = U(2:N+1);

%plots
figure(2);
plot(0:N, y_mpc, '-o', 'DisplayName', '$y_{mpc}$')
hold on;
plot(0:N-1, u_mpc, '-o', 'DisplayName', '$u_{mpc}$')
title('MPC system vs open loop system')
xlabel('timestep[n]')
legend('Interpreter', 'Latex', 'FontSize', 14);
