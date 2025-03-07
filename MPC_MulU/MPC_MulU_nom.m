%清屏
clc;
clear all;
close;

% 加载工具箱
% pkg load optim;

% 状态空间
% x(k+1) = A * x(k) + B * u(k)
% y(k) = C * x(k)

% 系统矩阵
A = [1 0.1; 0 2];
[n, ~] = size(A);
% 控制输入矩阵
B = [0 3;1 0.5];
[~, p] = size(B);
% 权重矩阵
Q = diag([0.1, 0.5]);
R = diag([1, 1]);
% 终端代价矩阵
F = diag([1, 1]);

% 设置step数量
steps = 100;
% 预测步长
N = 10;

% 初始化
x0 = [20; -20]
xk = zeros(n, steps);
uk = zeros(p, steps);
xk(:, 1) = x0;

% 损失函数中矩阵计算
[H, E, G, M, T] = MPC_MulU_CostMatrixSolve(A, B, Q, R, F, N);

% 循环
for k = 1:1:steps
    % 求解uk
    uk(:, k) = MPC_MulU_Prediction(xk(:, k), H, E, p);
    % 计算第k+1步的状态变量
    xk(:, k+1) = A * xk(:, k) + B * uk(:, k);
end

%% 绘图
subplot(2,1,1);
hold;
for i = 1:1:size(xk, 1)
    plot(xk(i, :));
end
legend("x");
hold off;

subplot(2, 1, 2);
hold;
for i = 1:1:size(uk, 1)
    plot(uk(i, :));
end
legend("u");




















