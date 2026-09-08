%% TPIM versus SA algorithm 
close all; 
%% 1. 输入待求的 N=100 的相互作用矩阵
N = 100;
J = Jmobius;
J = (J + J') / 2;

% 可选项，可以去掉对角线进行计算，也可不去(也即使用Cholesky分解后的)，默认带有对角线
% J = J - diag(diag(J)); 
%% 2. 算法梯度与步数设置
temstep = 80;              % 降温梯度数
step = 400;                 % 每个温度下的迭代步数
totstep = step * temstep;

%% ====================================================================
%% 3.1 3N 算法时域流，平均正向势垒标定温度
disp('正在基于 3N 算法标定初始温度 T0_3N...');

num_samples = 50;
mean_positive_dH_list = zeros(1, num_samples);

for s_idx = 1:num_samples
    s_init = 2 * randi([0, 1], 1, N) - 1;
    s_init_H = -s_init * J * s_init';
    
    % 1. 模拟 3N 算法的时域序列: S0 -> S_flip (0.15N) -> S0
    flip_indices = randperm(N, round(0.5 * N));
    flip_mask = false(1, N); flip_mask(flip_indices) = true;
    s_flip = s_init .* (1 - 2 * flip_mask);
    
    % 2. 评估单次 3N 时域流内全部 2N 个候选状态的能量增量
    dH_stream = zeros(1, 2 * N);
    for j = 1:N
        cand1 = [s_flip(1:j), s_init(j+1:N)];
        dH_stream(j) = (-cand1 * J * cand1') - s_init_H;
        
        cand2 = [s_init(1:j), s_flip(j+1:N)];
        dH_stream(N+j) = (-cand2 * J * cand2') - s_init_H;
    end
    
    % 3. 只对 3N 流内大于 0 的正向上升势垒求平均值
    positive_dH = dH_stream(dH_stream > 0);
    if ~isempty(positive_dH)
        mean_positive_dH_list(s_idx) = mean(positive_dH);
    end
end

% 4. T0_3N 定义为 3N 流内平均正向势垒对应的权重设置为0.5(为减小震荡可以设置更小的比如0.3)
T0_3N = -mean(mean_positive_dH_list)/log(0.5);
T0_3N = 16;  % 也可以认为给定合适的初温
%% 3.2 SA初始温度 T0 标定 
disp('正在进行SA标定初始温度 T0_single...');
num_samples = 50;
delta_H_single = [];

for idx = 1:num_samples
    s_temp = 2 * randi([0, 1], 1, N) - 1;

    % 单自旋算法抽样 (1个自旋翻转)
    flip_single = randi(N);
    s_probe_single = s_temp;
    s_probe_single(flip_single) = -s_temp(flip_single);
    dH_single = (-s_probe_single * J * s_probe_single') - (-s_temp * J * s_temp');
    if dH_single > 0
        delta_H_single = [delta_H_single, dH_single]; 
    end
end

% T0_single 定义为初始不良步的接受概率设为0.5(可以更小比如0.3)
T0_single = -mean(delta_H_single) / log(0.5);
% T0_single = 45; % 也可以人为设置合理的温度
%% 4：退火衰减率 dec 计算 (保证末期翻转自旋数量 k=1)
%% ====================================================================
dec = (2 / N)^(1 / (temstep - 1));
disp('--------------------------------------------------');
disp(['初始温度 T0 = ', num2str(T0_3N)]);
disp(['退火衰减率 dec = ', num2str(dec)]);
disp('--------------------------------------------------');

%% ====================================================================
%% 算法一：3N-SSflipS + 阶梯降温 + 类 Softmax 选取 (选择后接收率A = 1)
%% ====================================================================
disp('正在运行： 3N-SSflipS ...');
counum = 1;                              % 运行总次数，用于生成histo
Energy_method3N = zeros(counum,totstep); % 记录所有的迭代ising Energy H
for countnumber = 1:counum
    initial_sigma = 2 * randi([0, 1], 1, N) - 1;
    sigma_new = initial_sigma;
    H_current_new = -sigma_new * J * sigma_new'; % 随轮次改变的S0的H
    H_history_new = zeros(1, totstep);   % 在每一次内接收所有的H，最终取全局最小
    
    T = T0_3N;
    iter = 0;
    for circ = 1:temstep
        % 自旋转变数量 k (在当前 circ 迭代次数内为固定常数，随 T 降低而递减)
        k = max(1, round(0.5 * N * (T / T0_3N)));
        
        for p_idx = 1:step
            iter = iter + 1;
            
            % 1. 随机抽取 k 个自旋物理序号进行翻转
            flip_indices = randperm(N, k);
            flip_mask = false(1, N);
            flip_mask(flip_indices) = true;
            sigma_flip = sigma_new .* (1 - 2 * flip_mask);
            
            % 2. 产生 2N+1 个滑动窗口候选状态（包括初态S0）
            H_pool = zeros(1, 2 * N + 1);
            H_pool(1) = H_current_new; % 窗口 0 为当前状态 S0
            candidates = zeros(2 * N, N);

            % 前半段：S 逐渐过渡到 S_flip
            for j = 1:N
                candidates(j, :) = [sigma_flip(1:j), sigma_new(j+1:N)];
                H_pool(j+1) = -candidates(j, :) * J * candidates(j, :)';
            end
            % 后半段：S_flip 逐渐过渡回 S
            for k_idx = 1:N
                candidates(N+k_idx, :) = [sigma_new(1:k_idx), sigma_flip(k_idx+1:N)];
                H_pool(N+k_idx+1) = -candidates(N+k_idx, :) * J * candidates(N+k_idx, :)';
            end
            
            % Softmax 选择阶段：使用当前退火温度 T
            dH = H_pool - H_current_new;
            weights = exp(-dH / T); 
            weights(1) = 0; % 第一个为S0，不参与比较
            probs = weights / sum(weights); 
            
            % 投色子选择
            cum_probs = cumsum(probs);
            r = rand;
            sel_idx = find(cum_probs >= r, 1, 'first');
            
            % 接收阶段：A 恒等于 1，直接接受转移
            if sel_idx > 1
                sigma_new = candidates(sel_idx - 1, :);
                H_current_new = H_pool(sel_idx);
            end
            H_history_new(iter) = H_current_new;
        end
        % 温度阶梯衰减
        T = T * dec;
    end 
Energy_method3N(countnumber,:) = H_history_new;
end
Energy_minofmethod3N = min(Energy_method3N,[],2);
disp('3N-SSflipS + Softmax 算法运行完成！');
%% ====================================================================
%% 算法二：经典基准算法（标准单自旋翻转 Metropolis-Hastings）
% ====================================================================
disp('正在运行：SA单自旋翻转 退火...');
Tall = zeros(1,totstep);
counum = 1;
Energy_singlespin = zeros(counum,totstep);      % 记录所有的迭代ising Energy H
for countnumber = 1:counum
    initial_sigma = 2 * randi([0, 1], 1, N) - 1;
    sigma_classic = initial_sigma;
    H_current_classic = -sigma_classic * J * sigma_classic';
    H_history_classic = zeros(1, totstep);
    
    T = T0_single;
    iter = 0;
    for circ = 1:temstep
        for p_idx = 1:step
            iter = iter + 1;
            if countnumber == 1 
            Tall(iter) = T;
            end
            % 1. 随机选择一个自旋进行单点翻转
            flip_idx = randi([1, N]);
            sigma_probe = sigma_classic;
            sigma_probe(flip_idx) = -sigma_classic(flip_idx);
            
            % 2. 计算能量差
            H_probe = -sigma_probe * J * sigma_probe';
            dH = H_probe - H_current_classic;
            
            % 3. Metropolis 判定
            if dH < 0
                sigma_classic = sigma_probe;
                H_current_classic = H_probe;
            else
                p_accept = exp(-dH / T);
                if rand < p_accept
                    sigma_classic = sigma_probe;
                    H_current_classic = H_probe;
                end
            end       
            H_history_classic(iter) = H_current_classic;
        end
        T = T * dec;
    end
    Energy_singlespin(countnumber,:) = H_history_classic;
end
Energy_min_singlespin = min(Energy_singlespin,[],2);
disp('经典SA算法运行完成！');
%% ====================================================================
%% 5. 结果对比可视化
%% ====================================================================
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
plot(H_history_classic, 'Color', [0.85, 0.32, 0.09], 'LineWidth', 1.8, 'DisplayName', 'Classical SA');
hold on;
plot(H_history_new, 'Color', [0, 0.44, 0.74], 'LineWidth', 1.8, 'DisplayName', 'Proposed SSflipS (3N) with Softmax');
title(['Ising Ground State Search Comparison (N = ', num2str(N), ', Fully Adaptive)'], 'FontSize', 12);
xlabel('Monte Carlo Iteration Steps', 'FontSize', 11);
ylabel('Hamiltonian Energy', 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 10);
grid on;
hold off;
disp('--------------------------------------------------');
disp(['[经典 MH] 最终收敛能量: ', num2str(min(Energy_min_singlespin ))]);
disp(['[自适应 3N-SSflipS] 最终收敛能量: ', num2str(min(Energy_minofmethod3N))]);
disp('--------------------------------------------------');
%% 6. 绘制基态能量直方分布图(SA)
% 自动提取所有出现过的能量值和出现频次
[counts, centers] = groupcounts(Energy_min_singlespin); 
labels = string(centers); 
fig2 = figure('Color', 'w');
% 调整图形窗口大小
set(fig2, 'Units', 'centimeters', 'Position', [5, 5, 11.5, 8.5]);

% 绘制低饱和度蓝色半透明柱状图
bar_width = 0.35;
b1 = bar(categorical(labels, labels), counts, bar_width, ...
         'FaceColor', [0.12, 0.53, 0.90], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Hamiltonian Energy', 'FontSize', 10, 'FontName', 'Times New Rome');
ylabel('Number of Runs', 'FontSize', 10, 'FontName', 'Times New Rome');
title(['single Success Probability Distribution (N = ', num2str(N), ')'], ...
      'FontSize', 11, 'FontName', 'Times New Rome', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Rome', 'FontSize', 9);
grid on;
%% 7. 绘制基态能量直方分布图（3N+softmax）
% 自动提取所有出现过的能量值和出现频次
[counts, centers] = groupcounts(Energy_minofmethod3N); 
labels = string(centers); 

fig3 = figure('Color', 'w');
% 调整图形窗口大小
set(fig3, 'Units', 'centimeters', 'Position', [5, 5, 11.5, 8.5]);

% 绘制低饱和度蓝色半透明柱状图
bar_width = 0.35;
b2 = bar(categorical(labels, labels), counts, bar_width, ...
         'FaceColor', [0.12, 0.53, 0.90], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

xlabel('Hamiltonian Energy', 'FontSize', 10, 'FontName', 'Times New Rome');
ylabel('Number of Runs', 'FontSize', 10, 'FontName', 'Times New Rome');
title(['3N Success Probability Distribution (N = ', num2str(N), ')'], ...
      'FontSize', 11, 'FontName', 'Times New Rome', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Rome', 'FontSize', 9);
grid on;
%% 保存相关数据
% 保存算法对比
% save('ATPIMrand_runs100_data_Jmobius80400.mat',"dec","T0_single","T0_3N",'Jmobius', 'Energy_method3N','Energy_minofmethod3N','Energy_min_singlespin', 'Energy_singlespin');
% 保存纯SA结果用于和光学仿真对比
% save('ATPIMsingleflip_runs100_data_Jrank55040',"dec","T0_single","Tall",'Jrank5', 'Energy_min_singlespin', 'Energy_singlespin');