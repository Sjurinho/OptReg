%% Constants and system matrices
k1 = 1; k2 = 1; k3 = 1; T = 0.1;
A = [0 1;
    -k2 -k1];
B = [0 k3]';
C = [1 0];

x0 = [5 1]'; x0_hat = [6 0]';

%Discretized

Ad = eye(2) + A*T;
Bd = B*T;

%% Solve problem
Q = [4 0;
    0 4];
R = 1;
%% PROBLEM 3 - EXTENDED TO USE MPC
M = 10; %LQR timestep
N = 50; %timesteps
nx = 2; %state dimension
nu = 1; %input dimension
t = 0:T:N*T;

G_left = kron(eye(M), Q);
G_right = kron(eye(M), R);
G = blkdiag(G_left, G_right);

Aeq_left = kron(eye(M), eye(nx));
for i = 3:nx:M*nx
   Aeq_left(i:i+nx-1, i-nx:i-1) = -Ad;
end
Aeq_right = kron(eye(M), -Bd);
Aeq = [Aeq_left, Aeq_right];

beq = zeros(M*nx, 1);

xu = Inf*ones(M*nx,1);
xl = -Inf*ones(M*nx, 1);
uu = 4*ones(M*nu,1);
ul = -4*ones(M*nu,1);

ub = [xu;uu];
lb = [xl;ul];
%%Simulate

%Dynamics
x_tnext = @(x_t, u) Ad*x_t + Bd*u;
x_tnext_hat = @(x_t_hat, x_t, u) Ad*x_t_hat + Bd*u + K_f*(C*x_t - C*x_t_hat);

x = zeros((N+1)*nx, 1);
x_hat = zeros((N+1)*nx, 1);
u = zeros(N*nu, 1);

%Start conditions
x(1:nx) = x0;
x_hat(1:nx) = x0_hat;
j = 1;
for i = nx+1:nx:(N+1)*nx
    beq(1:nx) = Ad*x_hat(i-nx:i-1);

    %Constraints
    z = quadprog(G, [], [], [], Aeq, beq, lb, ub);
        
    x(i:i+nx-1) = x_tnext(x(i-nx:i-1), z(M*nx+1));
    x_hat(i:i+nx-1) = x_tnext_hat(x_hat(i-nx:i-1), x(i-nx:i-1), z(M*nx+1));
    u(j) = z(M*nx+1);
    j = j+1;
end
%Extract solution
x1 = x(1:nx:(N+1)*nx);
x2 = x(2:nx:(N+1)*nx);
x1_hat = x_hat(1:nx:(N+1)*nx);
x2_hat = x_hat(2:nx:(N+1)*nx);

%% Plots
figure(2)
plot(0:T:(N-1)*T, u)
figure(1)
plot(t, x1, 'r', 'DisplayName', '$x1$')
hold on; grid on;
plot(t, x1_hat, '--r', 'DisplayName', '$\hat{x1}$')
plot(t, x2, 'b', 'DisplayName', '$x2$')
plot(t, x2_hat, '--b', 'DisplayName', '$\hat{x2}$')
legend('Interpreter', 'Latex')
title('States vs state estimators', 'Interpreter', 'Latex')
xlabel('$t$[sec]', 'Interpreter', 'Latex')