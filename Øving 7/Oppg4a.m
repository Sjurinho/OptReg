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
%% PROBLEM 4a
[K, P, E] = dlqr(Ad, Bd, 1/2*Q, 1/2*R);
K_f = place(Ad', C', [0.5 + 0.03j, 0.5 - 0.03j]')';

%% Simulate
N = 50; %timesteps
nx = 2; %state dimension
nu = 1; %input dimension
t = 0:T:(N-1)*T;

%Dynamics
x_tnext = @(x_t, x_t_hat) Ad*x_t - Bd*K*x_t_hat;
x_tnext_hat = @(x_t_hat, x_t) Ad*x_t_hat - Bd*K*x_t_hat + K_f*(C*x_t - C*x_t_hat);

x = zeros(N*nx, 1);
x_hat = zeros(N*nx, 1);
u = zeros(N*nu, 1);

%Start conditions
x(1:nx) = x0;
x_hat(1:nx) = x0_hat;

for i = nx+1:nx:N*nx
    x(i:i+nx-1) = x_tnext(x(i-nx:i-1), x_hat(i-nx:i-1));
    x_hat(i:i+nx-1) = x_tnext_hat(x_hat(i-nx:i-1), x(i-nx:i-1));
end
%Extract solution
x1 = x(1:nx:N*nx);
x2 = x(2:nx:N*nx);
x1_hat = x_hat(1:nx:N*nx);
x2_hat = x_hat(2:nx:N*nx);

%% Plots
figure(1)
plot(t, x1, 'r', 'DisplayName', '$x1$')
hold on; grid on;
plot(t, x1_hat, '--r', 'DisplayName', '$\hat{x1}$')
plot(t, x2, 'b', 'DisplayName', '$x2$')
plot(t, x2_hat, '--b', 'DisplayName', '$\hat{x2}$')
legend('Interpreter', 'Latex')
title('States vs state estimators', 'Interpreter', 'Latex')
xlabel('$t$[sec]', 'Interpreter', 'Latex')