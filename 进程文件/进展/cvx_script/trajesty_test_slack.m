% CVX 代码实现 UAV 轨迹优化问题（满足新增条件）

% 仿真参数
M = 2; % 无人机数量
N = 180; % 时隙数目
Su = 100; % 无人机节点活动范围 (m)
Pk = 1.5; % 无人机节点发射功率 (W)
Pj = 1.5; % 干扰节点发射功率 (W)
sigma2 = 10^(-120/10); % 噪声功率 (W)
beta0 = 10^(-50/10); % 参考距离信道增益 (d0 = 1m)
q0 = [200, 800; 800, 800]; % 无人机初始位置 (m)
qj0 = [500, 500]; % 干扰节点初始位置 (m)
max_iter = 50; % 最大迭代次数
tol = 1e-3; % 收敛阈值

% 初始化变量
q = repmat(q0, [1, 1, N]); % 无人机位置 (M x 2 x N)
qj = zeros(1, 2, N); % 干扰节点位置 (1 x 2 x N)
qj(1, :, 1) = qj0; % 干扰节点初始位置
for n = 2:N
    qj(1, :, n) = qj(1, :, n-1) + [1, 0]; % 干扰源沿 x 方向以 1m/帧 移动
end
z_ik = zeros(M, M, N); % 松弛变量 z_ik[n] = ||q_i[n] - q_k[n]||^2
z_ij = zeros(M, N); % 松弛变量 z_ij[n] = ||q_i[n] - q_j[n]||^2
y_ik = zeros(M, M, N); % 松弛变量 y_ik[n] = Pk * h_ik[n] / (Pj * h_ij[n] + sigma2)

% 迭代优化
for iter = 1:max_iter
    % 固定 q_i[n], 优化 z_ik[n], z_ij[n], y_ik[n]
    cvx_begin
        variables z_ik(M, M, N) z_ij(M, N) y_ik(M, M, N)
        maximize sum(sum(sum(log(1 + y_ik) / log(2))))
        subject to
            % 信道增益约束
            for n = 1:N
                for i = 1:M
                    for k = 1:M
                        if i ~= k
                            % 一阶泰勒展开近似 h_ik[n]
                            h_ik_approx = beta0 / z_ik(i, k, n);
                            h_ik_approx_taylor = beta0 / z_ik(i, k, n) - beta0 / (z_ik(i, k, n)^2) * (z_ik(i, k, n) - z_ik(i, k, n));
                            % 一阶泰勒展开近似 h_ij[n]
                            h_ij_approx = beta0 / z_ij(i, n);
                            h_ij_approx_taylor = beta0 / z_ij(i, n) - beta0 / (z_ij(i, n)^2) * (z_ij(i, n) - z_ij(i, n));
                            % 约束条件
                            y_ik(i, k, n) <= (Pk * h_ik_approx_taylor) / (Pj * h_ij_approx_taylor + sigma2);
                        end
                    end
                    z_ij(i, n) >= norm(q(i, :, n) - qj(:, n))^2;
                end
            end
            % 无人机位置约束
            for n = 1:N
                for i = 1:M
                    norm(q(i, :, n) - q0(i, :)) <= Su;
                end
            end
    cvx_end

    % 固定 z_ik[n], z_ij[n], y_ik[n], 优化 q_i[n]
    cvx_begin
        variables q(M, 2, N)
        maximize sum(sum(sum(log(1 + y_ik) / log(2))))
        subject to
            % 信道增益约束
            for n = 1:N
                for i = 1:M
                    for k = 1:M
                        if i ~= k
                            z_ik(i, k, n) >= norm(q(i, :, n) - q(k, :, n))^2;
                        end
                    end
                    z_ij(i, n) >= norm(q(i, :, n) - qj(:, n))^2;
                end
            end
            % 无人机位置约束
            for n = 1:N
                for i = 1:M
                    norm(q(i, :, n) - q0(i, :)) <= Su;
                end
            end
            % 无人机节点前后位置差控制在 1m 以内
            for n = 2:N
                for i = 1:M
                    norm(q(i, :, n) - q(i, :, n-1)) <= 1;
                end
            end
    cvx_end

    % 判断收敛
    if iter > 1 && abs(cvx_optval - prev_optval) < tol
        break;
    end
    prev_optval = cvx_optval;
end


% 输出结果
disp('优化完成！');
disp('无人机最优轨迹：');
disp(q);
disp('总传输速率：');
disp(cvx_optval);