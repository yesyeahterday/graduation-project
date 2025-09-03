% 参数初始化
beta0 = 1e-5;       % 参考信道增益 (公式7中的β0,-50dB = 1e-5)
sigma2 = 1e-12;     % 噪声功率 (公式9中的σ²)
P_j = 4;            % 干扰功率 (公式11中的Pj)
S_u = 6 * 5;            % 初始位置最大偏移 (公式3-4中的Su)
N_m = 1 * 5;        %单步移动范围
N = 80;             % 时隙数
M = 2;              % 无人机数量-
C = 5;              % 信道数
max_iter = 200;      % 最大交替优化迭代次数
tol = 1e-5;         % 收敛阈值

% 初始化无人机位置为指定坐标
initial_positions = [30, 30; 
                    30, 70; 
                    70, 30; 
                    70, 70] * 5;  % M=4架无人机的初始坐标

% 设定每个节点 i 的初始功率（可以自定义不同的初始值）
p_initial = [1; 1; 1; 1];  % 维度 M × 1，每个节点的初始功率
p = repmat(p_initial,1, N);  

% 每个时隙的初始位置均设置为指定坐标
q = repmat(initial_positions, [1, 1, N]);  % 维度：M×2×N，所有时隙初始位置相同
q_prev = q;                                % 保存上一轮迭代的位置
disp(q);

% 干扰源位置（假设初始位置为 [10,10]，向右匀速移动，每步+1）
q_j = zeros(1, 2, N);  % 维度为 1×2×N
for n = 1:N
    q_j(1, :, n) = [10 + n - 1 , 50] * 5;  % 每个时隙向右移动 1 个单位
end


eps = 1e-6;

% 交替优化主循环
for iter = 1:max_iter
    total_rate = 0;
    for l = 1:M     % 对每架无人机进行交替优化
        % 固定其他无人机的位置
        q_fixed = q;
        q_fixed(l,:,:) = [];
        
        % 构建当前无人机l的优化问题
        cvx_clear;
        cvx_solver mosek;
        cvx_begin %quiet
            variable q_l(2, N)       % 无人机l的轨迹（x,y坐标）
            variable z_l(N)          % 引入的辅助变量z_l[n]（公式28）
            expressions obj(N)       % 目标函数
            expressions phi_norm(N)  % 二次项的一阶泰勒近似（公式26）
            
            for n = 1:N
                % -------------------------------
                % 公式23：泰勒展开R_u1的线性近似
                % -------------------------------
                % 计算参考点q_l^{(r)}[n]（上一轮迭代的位置）
                q_l_prev = squeeze(q_prev(l,:,n))';
                
                % 计算梯度项（公式23）
                grad_term = zeros(2,1);
                for k = 1:M
                    if k ~= l
                        % 计算梯度项中的分母D_lk（公式23中的D_lk）
                        d_lk_prev = norm(q_l_prev - squeeze(q_prev(k,:,n)));
                        d_lj_prev = norm(q_l_prev - squeeze(q_j(1,:,n)));
                        h_lk_prev = beta0 / (d_lk_prev^2 + eps); % 避免除以0
                        h_lj_prev = beta0 / (d_lj_prev^2 + eps);
                        D_lk = p(k,n) * h_lk_prev + P_j *h_lj_prev + sigma2;

                        delta_l_k = q_l_prev - squeeze(q_prev(k,:,n))';
                        delta_l_j = q_l_prev - squeeze(q_j(1,:,n))';
                        grad_term = grad_term + (-2*beta0 / log(2)) * ( (p(k,n) * delta_l_k) / (norm(delta_l_k)^4 * D_lk) + P_j * delta_l_j / (norm(delta_l_j)^4 * D_lk)  );
                    end
                end
                for i = 1:M
                    if i ~= l
                        % 计算梯度项中的分母D_il（公式24中的D_il）
                        d_il_prev = norm(q_l_prev - squeeze(q_prev(i,:,n)));
                        d_ij_prev = norm(squeeze(q_prev(i,:,n)) - q_j(1,:,n));
                        h_il_prev = beta0 / (d_il_prev^2 + eps);
                        h_ij_prev = beta0 / (d_ij_prev^2 + eps);
                        D_il = p(l,n) * h_il_prev + P_j *h_ij_prev + sigma2;

                        delta_l_i = q_l_prev - squeeze(q_prev(i,:,n))';
                        grad_term = grad_term + (-2 * beta0 / log(2)) * (p(l,n) * delta_l_i) / (norm(delta_l_i)^4 * D_il);
                    end
                end

                % 泰勒展开后的R_u1近似（公式23）
                R_u1_approx = grad_term' * (q_l(:,n) - q_l_prev);
                
                % -------------------------------
                % 公式26-28：处理R_u2的非凸项
                % -------------------------------
                % 二次项的一阶泰勒展开（公式26）
                delta_l_j = q_l_prev - q_j(1,:,n)';
                phi_norm(n) = norm(delta_l_j)^2 + 2 * delta_l_j' * (q_l(:,n) - q_l_prev);
                
                % 引入约束：e^{z_l[n]} >= P_j / (phi_norm * sigma2) （公式28）
                e_z = exp(z_l(n));
                e_z <= (phi_norm(n) * sigma2) / P_j;  % 修正分母项
                
                % -------------------------------
                % 总目标函数（公式30）
                % -------------------------------
                % 忽略与当前无人机l无关的项（i≠l时的求和）
                obj(n) = R_u1_approx - M * log( 1 / e_z + sigma2) / log(2);
            end
            
            maximize sum(obj)
            subject to
                % 移动步长约束（公式1-2）
                abs(q_l(1,1) - initial_positions(l,1)) <= N_m; 
                abs(q_l(2,1) - initial_positions(l,2)) <= N_m; 
                for n = 2:N
                    abs(q_l(1,n) - q_l(1,n-1)) <= N_m;
                    abs(q_l(2,n) - q_l(2,n-1)) <= N_m;
                end
                % 初始位置偏移约束（公式3-4）
                for n = 1:N
                    abs(q_l(1,n) - initial_positions(l,1)) <= S_u;
                    abs(q_l(2,n) - initial_positions(l,2)) <= S_u;
                end
        cvx_end

        disp(q_l(:,n));
        disp(q_l_prev);
        disp(q_l(:,n) - q_l_prev);
        
        % 更新当前无人机的位置
        q(l,:,:) = reshape(q_l, 1, 2, N);
        total_rate = total_rate + cvx_optval;
    end


    % 检查收敛条件
    if iter > 1 && abs(total_rate - prev_rate) < tol
        break;
    end
    prev_rate = total_rate;
    q_prev = q;
end

% 输出结果
disp(['算法收敛于迭代次数：', num2str(iter)]);
disp('最终无人机轨迹：');
disp(q);

generate_trajectory_gif(q, q_j, initial_positions, S_u, N_m);

function generate_trajectory_gif(q, q_j, initial_positions, S_u, N_m)
    % 坐标缩放（原始数据需除以5）
    scale_factor = 5;
    q = q / scale_factor;
    q_j = q_j / scale_factor;
    initial_positions = initial_positions / scale_factor;
    S_u = S_u / scale_factor;
    N_m = N_m / scale_factor;
    
    % 创建画布
    fig = figure('Position', [100 100 800 600]);
    hold on;
    grid on;
    axis equal;
    
    % 坐标范围
    xlim([0, 100]);
    ylim([0, 100]);
    
    % 颜色定义
    colors = lines(size(q,1));
    
    % 绘制静态元素（移动范围正方形）
    for m = 1:size(q,1)
        x_center = initial_positions(m,1);
        y_center = initial_positions(m,2);
        rectangle('Position', [x_center-S_u, y_center-S_u, 2*S_u, 2*S_u],...
                 'EdgeColor', [colors(m,:) 0.3], 'LineStyle', '--', 'LineWidth', 1);
    end
    
    % 初始化动态元素
    h_drones = gobjects(size(q,1), 1); % 无人机位置标记句柄
    h_jammer = plot(NaN, NaN, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5); % 干扰源标记
    
    % 初始化无人机标记
    for m = 1:size(q,1)
        h_drones(m) = plot(NaN, NaN, 'o', 'Color', colors(m,:),...
            'MarkerSize', 8, 'MarkerFaceColor', colors(m,:));
    end
    
    % 设置标注
    xlabel('X 坐标 (米)');
    ylabel('Y 坐标 (米)');
    title('无人机及干扰源实时位置');
    %legend([h_drones; h_jammer], [arrayfun(@(m)sprintf('无人机%d',m)), 1:size(q,1), '干扰源']);
    
    % 生成GIF
    gif_filename = 'drone_trajectory.gif';
    for n = 1:size(q,3)
        % 更新无人机位置
        for m = 1:size(q,1)
            set(h_drones(m), 'XData', q(m,1,n), 'YData', q(m,2,n));
        end
        
        % 更新干扰源位置
        set(h_jammer, 'XData', q_j(1,1,n), 'YData', q_j(1,2,n));
        
        % 捕获帧
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im, 256);
        
        % 写入GIF
        if n == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
    
    close(fig);
    disp(['动态图已保存为: ', gif_filename]);
end



