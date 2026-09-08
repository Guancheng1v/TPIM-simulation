%% 1. 数据源
J = Jrank5;
J_num = size(J, 1);
if mod(J_num, 2) == 0
    J1gpu = zeros(J_num+1, J_num+1, 'gpuArray');
    J1gpu(1:end-1, 1:end-1) = gpuArray(J);
else
    J1gpu = gpuArray(J);
end
Hdiag = trace(J);
n = size(J1gpu, 1); % 自旋规模
disp(['系统自旋物理规模 N = ', num2str(n)]);

%% 2. Cholesky 矩阵分解，获取权重矩阵 Ugpu

% Cholesky分解，输出Xi矩阵和实际考虑的相互作用矩阵Jtest(加了对角线)
[U,Jtest] = JtoKsi(n,J1gpu);
R = rank(U);
column_sums = sum(abs(U)); 

% 使用微小阈值（如 1e-12）来判断哪些列是有值的
is_nonzero = (column_sums > 1e-12); 
% 核心：只提取有值的列,将 U 矩阵的尺寸从 101x101 压缩为 101xR
U = U(:, is_nonzero);
Ugpu = gpuArray(U);
disp(['Cholesky 分解完成，系统解耦维度（秩）为: ', num2str(R)]);

%% 3. 初始退火温度 T0 计算 ===
disp('正在估算合理的初始退火温度 T0...');
num_samples = n;
delta_H_list = [];
for s_idx = 1:num_samples
    s_temp = 2 * randi([0, 1], 1, n, 'gpuArray') - 1;
    flip_mask = (rand(1, n, 'gpuArray') < 0.2); 
    s_probe = s_temp .* (1 - 2 * flip_mask);
    dH = (-s_probe * Jtest * s_probe') - (-s_temp * Jtest * s_temp');
    if dH > 0
        delta_H_list = [delta_H_list, dH];
    end
end
dH_avg = mean(gather(delta_H_list));
T0_estimated = -dH_avg / log(0.4); % 初始接受概率设为 40%
disp(['-> 自动标定的推荐初始温度 T0 = ', num2str(T0_estimated)]);
%% 4. 退火算法参数初始化
step = 30;                                % 每个温度梯度的迭代次数               
temstep = 40;                             % 降温梯度数量
dec = (2 / n)^(1 / (temstep  - 1));
disp(['退火衰减率 dec = ', num2str(dec)]);
totstep = step * temstep;      
T0 = 40;                                  % 为方便和算法比较，可以人为选取推荐温度附近的值
counum = 1;                                % 独立实验总次数（100次用于标准测试）

% 定义高维记录矩阵
Energy_min = zeros(counum, 1);
Energy_iteration = zeros(counum, totstep);
Tall = zeros(1,totstep);
% 定义实时 H-I 对齐记录数据（仅在单次实验内使用）
H_step = zeros(1, totstep);
I_step = zeros(1, totstep);

%% 5. 物理传播参数初始化
nt = 2^18;                                    % 时域采样点数
deltaf = 98e9;                                % 频梳谱线间隔
T_p = 1/deltaf;
dt = T_p/128;
t_tol = dt * (nt-1);
pulsenum = 1;                                 % 单个自旋调制对应的脉冲数
lambda = 1560e-9;                             % 中心波长
D = 2311;                                     % CFBG色散系数

TPIM = TPIM_simulation(nt, t_tol, deltaf, lambda, pulsenum, D);
TPIM.M = 1;
dt = t_tol / (nt - 1); 

%% 6. 时间原点取值校准
calib_sequence = zeros(1, 3*n, 'gpuArray');
calib_sequence(1 : n) = ones(1, n, 'gpuArray'); 
calib_xi = ones(1, n, 'gpuArray'); 
calib_init = calib_sequence .* [calib_xi, calib_xi, calib_xi];

TPIM.create_combs(n, ones(1, n, 'gpuArray'));
TPIM.create_modu_ampandphase(3 * n, calib_init, 0, 0);
TPIM.create_inputsignal();
TPIM.create_propagation(TPIM.D, 1);
[~, idx_win0] = max(abs(TPIM.Elt_new).^2);

%% 7. 自动计算数学哈密顿量与物理光强的转换系数 I2E
disp('正在自动标定转换系数 (I2E)...');
test_sigma = ones(1, n, 'gpuArray'); % 规定全一用于计算比值，如果为0可以人为选取其它状态
H_math_test = -test_sigma*Jtest*test_sigma';  
% 调用 3N 单次物理传播函数
I_windows_test = get_intensities_3N_single_pass(test_sigma, test_sigma, Ugpu, R, TPIM, deltaf, dt, nt, idx_win0);
H_opt_test = I_windows_test(1); 

I2E = H_opt_test / H_math_test;              
disp(['标定完成，I2E = ', num2str(I2E)]);

%% 8. 多次独立退火实验主循环
disp('==================================================');
disp(['开始运行 ', num2str(counum), ' 次独立全光子 Ising 退火实验...']);
tic; 

for count = 1:counum
    sigma = 2 * randi([0, 1], 1, n, 'gpuArray') - 1;   
    T = T0;
    circ = 0;                           
    Hamiltonmin = zeros(1, totstep, 'gpuArray');     
    
    % 实时控制台看板
    reverseStr = repmat('\b', 1, 0);
    
    for iter = 1:totstep
        % 1. 在当前温度阶梯下，确定翻转自旋个数 k
        k = max(1, round(0.5 * n * (T / T0)));
        
        % 计算内部步数 pin 
        pin = mod(iter-1, step) + 1;
        
        % 2. 随机无重复地抽取 k 个自旋序号进行翻转
        flip_indices = randperm(n, k);
        flip_mask = false(1, n);
        flip_mask(flip_indices) = true;
        sigma_flip = sigma .* (1 - 2 * gpuArray(flip_mask));
       
        % 单次 3N (SSflipS) 物理传播与时序采样
        I_windows = get_intensities_3N_single_pass(sigma, sigma_flip, Ugpu, R, TPIM, deltaf, dt, nt, idx_win0);
        H0 = I_windows(1); 

        %  Softmax 选择
        dH = I_windows - H0;
        weights = exp(-dH / (T * I2E)); 
        weights(1) = 0;
        probs = weights / sum(weights); 
        
        probs_cpu = gather(probs); 
        cum_probs = cumsum(probs_cpu); 
        r = rand; 
        sel_idx = find(cum_probs >= r, 1, 'first'); 
        
        % 记录当前比较选中的候选状态 H 与 I（仅第1次记录，用于 H-I 曲线）
        if count == 1 
            Tall(iter) = T;
        end
            if sel_idx > 1
                j_sel = sel_idx - 1;
                if j_sel <= n
                    sigma_candidate = [sigma_flip(1 : j_sel), sigma(j_sel + 1 : n)];
                else
                    k_sel = j_sel - n;
                    sigma_candidate = [sigma(1 : k_sel), sigma_flip(k_sel + 1 : n)];
                end
                if count == 1
                I_step(iter) = I_windows(sel_idx); 
                H_step(iter) = -sigma_candidate * Jtest * sigma_candidate'; 
                end
            else
                if count == 1
                I_step(iter) = H0;
                H_step(iter) = -sigma * Jtest * sigma';
                end
            end

        % 接收判定阶段（A = 1，直接接受转移）
        if sel_idx > 1
            sigma = sigma_candidate;
            H0 = I_windows(sel_idx); 
        end
        
        % 记录当前的数值计算 Hamiltonian
        Hamiltonmin(iter) = -sigma * Jtest * sigma';
        
        % 进度看板
        msg = sprintf('[实验 %d/%d] 迭代: %5d/%d',count, counum, iter, totstep);
        fprintf([reverseStr, msg]);
        reverseStr = repmat('\b', 1, length(msg));
      
        % 温度衰减（每个阶梯步数末尾执行）
        if pin == step
            circ = circ + 1;
            T = T * dec;
        end
    end
    
    % 数据记录
    Hamiltonmin_cpu = gather(Hamiltonmin);
    Energy_min(count, 1) = min(Hamiltonmin_cpu);         
    Energy_iteration(count, :) = Hamiltonmin_cpu;        
    fprintf('\n');
end

time_elapsed = toc; 
disp('==================================================');
disp(['独立实验全部结束！100次物理传播迭代计算总耗时: ', num2str(time_elapsed), ' 秒。']);
%% 9. 绘制 3N 算法实时 Hamiltonian 退火下降曲线（只展示第1次实验轨迹作为代表）

% 可选择使用不带有对角线的H，默认为带有对角线的
% Energy_iteration = Energy_iteration + Hdiag;
% H_step = H_step + Hdiag;
% Energy_min = Energy_min + Hdiag;
figure(1);
plot(Energy_iteration(1, :)', 'LineWidth', 1.2);
title(['Proposed SSflipS Hamiltonian Trajectories (N = ', num2str(n), ')']);
xlabel('Monte Carlo Iteration Steps');
ylabel('Hamiltonian Energy');
grid on;

%% 10. 绘制 H-I 散点图
figure(2);
scatter(H_step, I_step / I2E, 35, 'ob', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
% 绘制理想状态 y = x 的参考对角虚线
plot([min(H_step) max(H_step)], [min(H_step) max(H_step)] , 'r--', 'LineWidth', 2);
xlabel('H/J', 'FontName', 'Times New Rome', 'FontSize', 16);
ylabel('I/J', 'FontName', 'Times New Rome', 'FontSize', 16);
set(gca, 'FontSize', 16); % 设置坐标轴刻度字体
grid off;
box on;
hold off;

%% 11. 绘制基态能量直方图
% 提取所有出现过的能量值和出现频次
[counts, centers] = groupcounts(Energy_min); 
labels = string(centers); 

fig3 = figure('Color', 'w');
% 调整图形窗口大小
set(fig3, 'Units', 'centimeters', 'Position', [5, 5, 11.5, 8.5]);

% 绘图
bar_width = 0.35;
b1 = bar(categorical(labels, labels), counts, bar_width, ...
         'FaceColor', [0.12, 0.53, 0.90], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

xlabel('Hamiltonian Energy', 'FontSize', 10, 'FontName', 'Times New Rome');
ylabel('Number of Runs', 'FontSize', 10, 'FontName', 'Times New Rome');
title(['Success Probability Distribution (N = ', num2str(n), ')'], ...
      'FontSize', 11, 'FontName', 'Times New Rome', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Rome', 'FontSize', 9);
grid on;
%% 

% 保存当前统计数据
% save('TPIM3N_multi_data_Jrank55040.mat',"dec","T0","Tall",'Jtest','Jrank5', 'Energy_min', 'Energy_iteration',"I2E","H_step","I_step");
% disp('统计退火数据已成功保存');

%% 单次 3N 物理传播与时序对齐采样
function I_windows = get_intensities_3N_single_pass(sigma, sigma_flip, Ugpu, R, TPIM, deltaf, dt, nt, idx_win0)
    n = size(Ugpu, 1);
    % 3N 滑动窗口会产生 2N+1 个干涉光强采样点
    I_windows = zeros(1, 2*n+1, 'gpuArray');
    I_tot_waveform = zeros(1, nt, 'gpuArray');
    
    % 在时域拼接成 SSflipS (3N 长度) 的物理自旋状态
    total_sequence = [sigma, sigma_flip, sigma]; 
    
    % 产生光梳
    freqa = ones(1, n, 'gpuArray');
    TPIM.create_combs(n, freqa);
    
    % 一次性计算 R 个通道的 3N 物理调制与传播
    for rJ = 1:R
        xi_k = Ugpu(:, rJ)';                  
        xi_k_replicated = [xi_k, xi_k, xi_k]; % 进行 3N 长度复制
        
        % 调制信号
        initheav = total_sequence .* xi_k_replicated; 
        
        % 调制并传播
        % 调制器会自动在时域脉冲加载相应的交替相位补偿信号 (-1)^K
        TPIM.create_modu_ampandphase(3 * n, initheav, 0, 0);
        TPIM.create_inputsignal();
        TPIM.create_propagation(TPIM.D, 1);
        
        % 叠加非相干功率
        I_tot_waveform = I_tot_waveform + abs(TPIM.Elt_new).^2;
    end
    
    % =================================================================
    % 基于 Window 0 进行时序采样：
    % 采样间隔为 1/deltaf，总共采样 2N+1 个点
    % =================================================================
    for j = 0:(2*n)
        time_offset = j * (1 / deltaf); 
        idx_j = idx_win0 + round(time_offset / dt); 
        
        % 记录对应滑动窗口的光强
        I_windows(j+1) = -I_tot_waveform(idx_j); 
    end
end