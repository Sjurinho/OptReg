%%Model
A = [0 0 0;
    0 0 1;
    0.1 -0.79 1.78];
B = [1 0 0.1]';
C = [0 0 1];

%Objectives
Q = [0 0 0;
    0 0 0;
    0 0 2];
R = 2;

%Start condition
x0 = [0 0 1]';


%%Generate over time horizon N
N = 30;
nx = size(A,2);
mu = size(B,2);

%Bounds 
x_bounds = [Inf(N*nx,1) -Inf(N*nx,1)];
u_bounds = [ones(N*mu,1) -1*ones(N*mu,1)];
ub = [x_bounds(:,1); u_bounds(:,1)];
lb = [x_bounds(:,2); u_bounds(:,2)];

%Equality constraints
Aeq = gen_aeq(A, B, N, nx, mu); %From lab exercise

beq = zeros(1, N*nx);
beq(1,1:nx) = A*x0;

%Transform system into z vector
G = zeros(N*nx+N*mu, N*nx+N*mu);
for i = 1:nx:N*nx+1
    G(i:i+nx-1,i:i+nx-1) = blkdiag(Q);
end
for j = N*nx+1:N*nx + N*mu
    G(j, j) = blkdiag(R);
end
