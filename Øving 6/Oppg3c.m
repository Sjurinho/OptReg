%%Model
A = [0 0 0;
    0 0 1;
    0.1 -0.79 1.78];
B = [1 0 0.1]';
C = [0 0 1];

B_length = [1, 1, 2, 4, 8, 14]; %Lengths of input
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
u_bounds = [ones(N*mu/5,1) -1*ones(N*mu/5,1)];
ub = [x_bounds(:,1); u_bounds(:,1)];
lb = [x_bounds(:,2); u_bounds(:,2)];

%Equality constraints
ones_block = blkdiag(ones(B_length(1),1), ...
                     ones(B_length(2),1), ...
                     ones(B_length(3),1), ...
                     ones(B_length(4),1), ...
                     ones(B_length(5),1), ...
                     ones(B_length(6),1));      % Block-diagonal matrix of 1-vectors
A_input = kron(ones_block, -B);                  % Component 3 of A_eq
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

% Cost function
Qt = 2*diag([0, 0, 1]);
Q = kron(eye(N), Qt);
Rt = 2*1;
R = kron(diag(B_length), Rt); % Note that each R in G is multiplied by the block length!
G = blkdiag(Q, R);

[z, fval, exitflag, opts] = quadprog(G,[], [], [], Aeq, beq, lb, ub, x0); %1f
%Extract solution

y = [x0(3) z(nx:3:N*nx)']';
u_block = z(N*nx+1:N*nx+N/5*mu);
u = [];
for j=1:length(B_length)
    for i = 1:B_length(j)
        u = [u u_block(j)];
    end     
end
disp(opts);

%% Plots
figure(3)
plot(0:N, y, '-o')
hold on;
plot(0:N-1, u, '-ro')
title('Open loop optimized system')
xlabel('timestep[n]')
legend({'$y$', '$u$'}, 'Interpreter', 'Latex', 'FontSize', 14);
hold off;