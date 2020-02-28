%% Constants and system matrices
m = 1;
A = [1 0.5;
    0 1];
B = [0.125/m 0.5/m]';

%Control matrices - dlqr
Q = [2 0;
    0 2];
R = 2;

%% Solve Ricatti equation
[K, P, E] = dlqr(A, B, 1/2*Q, 1/2*R, []);