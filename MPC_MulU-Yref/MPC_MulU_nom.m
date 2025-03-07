%清屏
clc;
clear all;
close all;

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
C = [1 0; 0 1];
% 权重矩阵
Q = diag([1, 2]);
R = diag([1, 5]);
% 终端代价矩阵
F = diag([5, 1]);

% 设置step数量
steps = 500;
% 预测步长
N_step = 20;

% 设置yr期望值
y_ref = [10; -5];  % 期望横摇角0° 期望航向角30°

% 状态空间 x(k+1) = A * x(k) + B * u(k)；  y(k) = C * x(k)
[nA, ~] = size(A);
[~, pB] = size(B);

% 初始化
x0 = [20; -20];
xk = zeros(nA, steps);
uk = zeros(pB, steps);
xk(:, 1) = x0;

% 损失函数中矩阵计算
[R_bar, Q_bar, PHI, Gamma, M, T] = MPC_MulU_CostMatrixSolve(A, B, C, Q, R, F, N_step);

% 循环
for k = 1:1:steps
    % 构造参考轨迹在预测时域内
    Y_ref = repmat(y_ref, N_step+1, 1); % 若是与时间相关的函数，修改此处
    % 求解uk
    uk(:, k) = MPC_MulU_Prediction(xk(:, k), pB, Y_ref, R_bar, Q_bar, PHI, Gamma);
    % 计算第k+1步的状态变量
    xk(:, k+1) = A * xk(:, k) + B * uk(:, k) + rand(nA, 1); % 添加扰动
end

%% 绘图
figure(1)
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



















