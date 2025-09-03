clc; clear; close all;

% 初始化参数
N = 180; % 时隙数目
M = 2; % 无人机数量
M_scope = 50; %无人机运动范围
M_distance = 1; %单步距离
q_initial = [200, 800; 800, 800]'; % 无人机初始位置 (2xM)
q_jammer_start = [500; 500]; % 干扰节点初始位置
P_k = 1.5; % 无人机发射功率
P_j = 1.5; % 干扰节点发射功率
sigma2 = 10^(-12); % 噪声功率
beta0 = 10^(-5); % 信道增益常数
max_iter = 100; % 最大迭代次数
tol = 1e-12; % 收敛阈值

% 初始化无人机位置
q = repmat(q_initial, [1, 1, N]); % q(i, m, n) -> 第n个时隙，无人机m的位置 (2D)

% 初始化 SCA 迭代
q_old = q;

for iter = 1:max_iter
    % 计算当前的信道增益
    h_ik = zeros(M, M, N); % 无人机之间的信道增益
    h_ij = zeros(M, N); % 无人机到干扰节点的信道增益

    for n = 1:N
        q_jammer = q_jammer_start + [n; 0]; % 干扰节点沿x轴每步移动1m

        for i = 1:M
            for k = 1:M
                if k ~= i
                    h_ik(i, k, n) = beta0 / norm(q(:, i, n) - q(:, k, n))^2;
                end
            end
            h_ij(i, n) = beta0 / norm(q(:, i, n) - q_jammer)^2;
        end
    end

    % 使用CVX求解优化问题
    cvx_begin quiet
        variable q_new(2, M, N) % 新的无人机位置
        
        % 目标函数：最大化信道容量
        expr = 0;
        for n = 1:N
            for i = 1:M
                interference = P_j * h_ij(i, n) + sigma2; % 干扰加噪声

                expr_i = sum(log (1 + (P_k * h_ik(i, :, n)) / interference));
                expr = expr + expr_i;
            end
        end
        
        maximize expr / log(2);
        
        subject to
            % 约束 1: 无人机不能超出初始位置 50m
            for n = 1:N
                for i = 1:M
                    norm(q_new(:, i, n) - q_initial(:, i)) <= M_scope;
                end
            end
            
            % 约束 2: 无人机每一步的移动距离不能超过 1m
            for n = 2:N
                for i = 1:M
                    norm(q_new(:, i, n) - q_new(:, i, n-1)) <= M_distance;
                end
            end
    cvx_end

    % 检查收敛条件
    if norm(q_new - q_old, 'fro') < tol
        disp(['Converged at iteration ', num2str(iter)]);
        break;
    end
    
    % 更新无人机位置
    q_old = q;
    q = q_new;
end

% 绘制无人机和干扰节点的轨迹
figure; hold on; grid on;
colors = ['r', 'b'];
for i = 1:M
    plot(squeeze(q(1, i, :)), squeeze(q(2, i, :)), [colors(i), '-o'], 'LineWidth', 1.5);
end
plot(500:500+N-1, repmat(500, 1, N), 'k--', 'LineWidth', 2); % 干扰节点轨迹
legend('无人机 1', '无人机 2', '干扰节点');
xlabel('X 坐标 (m)');
ylabel('Y 坐标 (m)');
title('无人机轨迹优化');
