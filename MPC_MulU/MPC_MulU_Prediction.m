%% 预测求解uk
function uk = MPC_MulU_Prediction(xk0, H, E, p)
    f = xk0' * E;
    U_k = quadprog(H, f); % E的矩阵大小固定住了预测步长  %% 可以添加约束
    uk = U_k(1:p, 1); % 取第一个预测结果
end
