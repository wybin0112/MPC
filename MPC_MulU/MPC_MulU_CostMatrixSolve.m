%% 损失函数中的矩阵计算
function [H, E, G, M, T] = MPC_MulU_CostMatrixSolve(A, B, Q, R, F, N)

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

% 计算新权值矩阵Q_bar
% [nQ, mQ] = size(Q);
% [nF, mF] = size(F);
% Q_bar = zeros(N*nQ+nF, N*mQ+mF);
% for i = 0:1:N-1
%     Q_bar((i*nQ+1):(i+1)*nQ, (i*nQ+1):(i+1)*nQ) = Q;
% end
% Q_bar((N*nQ+1):end, (N*nQ+1):end) = F;
Q_bar = kron(eye(N), Q);
Q_bar = blkdiag(Q_bar, F);


% 计算新权值矩阵R_bar
% [nR, mR] = size(R);
% R_bar = zeros(N*nR, N*mR);
% for i = 0:1:N-1
%     R_bar((i*nR+1):(i+1)*nR, (i*nR+1):(i+1)*nR) = R;
% end
R_bar = kron(eye(N), R);

% 矩阵计算结果
H = (T') * Q_bar * T + R_bar;
E = (M') * Q_bar * T;
G = (M') * Q_bar * M;

end
