%% Constants and system matrices
k1 = 1; k2 = 1; k3 = 1; T = 0.1;
A = [0 1;
    -k2 -k1];
B = [0 k3]';
C = [1 1];

x0 = [5 1]';

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

G_left = kron(eye(M-1), Q);
G_right = kron(eye(M), R);
G = blkdiag(G_left,2*P, G_right);

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

x = zeros((N+1)*nx, 1);
u = zeros(N*nu, 1);
x_ol = zeros((N+1)*nx,1);

%Start conditions
x(1:nx) = x0;
j = 1;
for i = nx+1:nx:(N+1)*nx
    beq(1:nx) = Ad*x(i-nx:i-1);

    %Constraints
    z = quadprog(G, [], [], [], Aeq, beq, lb, ub);
        
    x(i:i+nx-1) = x_tnext(x(i-nx:i-1), z(M*nx+1));
    u(j) = z(M*nx+1);
    j = j+1;
end
%Extract solution
x1 = x(1:nx:(N+1)*nx);
x2 = x(2:nx:(N+1)*nx);

%% Plots
figure(2)
plot(0:T:(N-1)*T, u) 
grid on;
figure(1)
plot(t, x1, 'r', 'DisplayName', '$x1$')
hold on; grid on;
plot(t, x2, 'b', 'DisplayName', '$x2$')
legend('Interpreter', 'Latex')
title('States vs state estimators', 'Interpreter', 'Latex')
xlabel('$t$[sec]', 'Interpreter', 'Latex')