%% 损失函数中的矩阵计算
function [R_bar, Q_bar, PHI, Gamma, M, T] = MPC_MulU_CostMatrixSolve(A, B, C, Q, R, F, N)

% X_0~N = M * x0 + T * U_0~(N-1)

% 计算M矩阵
[nA, ~] = size(A);
M = zeros(nA*(N+1), nA);
for i = 0:1:N
    A_temp = A^i;
    M((i*nA+1):(i+1)*nA, 1:nA) = A_temp;
end

% 计算T矩阵
[nB, mB] = size(B);
T = zeros(nB*(N+1), mB*N);
for i = 0:1:N-1
    for j = i:-1:0
        AjB = (A^j) * B;
        T(((i+1)*nB+1):(i+2)*nB, ((i-j)*mB+1):(i-j+1)*mB) = AjB;
    end
end

% 计算新权值矩阵Q_bar--Y
Q_bar = kron(eye(N), Q);
Q_bar = blkdiag(Q_bar, F);

% 计算新权值矩阵R_bar--U
R_bar = kron(eye(N), R);

% Y = PHI*xk + Gamma*U;
[nC, mC] = size(C);
PHI = zeros(nC*(N+1), mC);
for i = 0:1:N
    C_temp = C*A^i;
    PHI((i*nC+1):(i+1)*nC, 1:mC) = C_temp;
end

Gamma = zeros(nC*(N+1), mB*N);
for i = 0:1:N-1
    for j = i:-1:0
        CAjB = C*(A^j) * B;
        Gamma(((i+1)*nC+1):(i+2)*nC, ((i-j)*mB+1):(i-j+1)*mB) = CAjB;
    end
end

end
