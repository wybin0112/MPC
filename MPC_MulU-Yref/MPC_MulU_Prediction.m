%% 预测求解uk
function uk = MPC_MulU_Prediction(xk0, p, Yr, R_bar, Q_bar, PHI, Gamma)
    H = Gamma' * Q_bar * Gamma + R_bar;
    f = Gamma' * Q_bar * (PHI * xk0 - Yr); % 列向量

    % 输入约束
    lb = [];
    ub = [];
    U_k = quadprog(H, f, [],[],[],[],lb,ub);

    uk = U_k(1:p, 1); % 取第一个预测结果
end

