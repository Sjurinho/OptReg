%%Model
A = [0 0 0;
    0 0 1;
    0.1 -0.79 1.78];
B = [1 0 0.1]';
C = [0 0 1];

B_length = 5;%Block length

%Objectives
Q = [0 0 0;
    0 0 0;
    0 0 2];
R = 2;

%Start condition
x0 = [0 0 1]';


%% Generate over time horizon N
N = 30;
nx = size(A,2);
mu = size(B,2);

%Bounds 
x_bounds = [Inf(N*nx,1) -Inf(N*nx,1)];
u_bounds = [ones(N*mu/B_length,1) -1*ones(N*mu/B_length,1)];
ub = [x_bounds(:,1); u_bounds(:,1)];
lb = [x_bounds(:,2); u_bounds(:,2)];

%Equality constraints
A_input = kron(kron(eye(6), ones(5,1)),-B); %The input part
A_state = zeros(N*nx); %the state part
A_state(1:1+nx-1, 1:1+nx-1) = eye(nx);
for i = nx+1:nx:N*nx
    A_state(i:i+nx-1, i:i+nx-1) = eye(nx);
    A_state(i:i+nx-1, i-nx:i-1) = -A;
end
%Combine
Aeq = [A_state A_input];

beq = zeros(1, N*nx);
beq(1,1:nx) = A*x0;

%Transform system into z vector
G = zeros(N*nx+N*mu/B_length, N*nx+N*mu/B_length);
for i = 1:nx:N*nx+1
    G(i:i+nx-1,i:i+nx-1) = blkdiag(Q);
end
for j = N*nx+1:N*nx + (N*mu/B_length)
    G(j, j) = B_length*blkdiag(R);
end

%% Plots
figure(3)
plot(0:N, y, '-o')
hold on;
plot(0:N-1, u, '-ro')
title('Open loop optimized system')
xlabel('timestep[n]')
legend({'$y$', '$u$'}, 'Interpreter', 'Latex', 'FontSize', 14);
hold off;