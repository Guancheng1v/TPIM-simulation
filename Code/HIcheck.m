%% 1. 初始化矩阵与物理参数
J = J22;
J_num = size(J, 1);
if mod(J_num, 2) == 0
    J1gpu = zeros(J_num+1, J_num+1, 'gpuArray');
    J1gpu(1:end-1, 1:end-1) = gpuArray(J);
else
    J1gpu = gpuArray(J);
end
n = size(J1gpu, 1); 

% Cholesky 分解
[Ksi, Jtest] = JtoKsi(n, J1gpu);
R = rank(Ksi); 

column_sums = sum(abs(Ksi)); 
Hdiag = trace(Jtest);
% 使用一个微小的阈值（如 1e-12）来判断哪些列是有值的
is_nonzero = (column_sums > 1e-12); 

% 核心：只提取有值的列。这会将 U 矩阵的尺寸从 101x101 压缩为 101xR！
Ksi = Ksi(:, is_nonzero);
Ugpu = gpuArray(Ksi);
%% 2. 色散传播
nt = 2^18;                                    
deltaf = 98e9;                                
T_p = 1/deltaf;
dt = T_p/128;
t_tol = dt * (nt-1);
pulsenum = 1;                                 
lambda = 1560e-9;                             
D = 2311;                                     

TPIM = TPIM_simulation(nt, t_tol, deltaf, lambda, pulsenum, D);
TPIM.M = 1;
dt = t_tol / (nt - 1); 

% 标定窗口取点 Window 0
calib_sequence = zeros(1, 3*n, 'gpuArray');
calib_sequence(1 : n) = ones(1, n, 'gpuArray'); 
calib_xi = ones(1, n, 'gpuArray'); 
calib_init = calib_sequence .* [calib_xi, calib_xi, calib_xi];
TPIM.create_combs(n, ones(1, n, 'gpuArray'));
TPIM.create_modu_ampandphase(3 * n, calib_init, 0, 0);
TPIM.create_inputsignal();
TPIM.create_propagation(TPIM.D, 1);
[~, idx_win0] = max(abs(TPIM.Elt_new).^2);

% 计算标定系数 I2E
test_sigma = ones(1, n, 'gpuArray');          
H_math_test = -test_sigma*Jtest*test_sigma';  
I_windows_test = get_intensities_3N_single_pass(test_sigma, test_sigma, Ugpu, R, TPIM, deltaf, dt, nt, idx_win0);
I2E = I_windows_test(1) / H_math_test;              

%% 3. 随机输入不同的自旋状态，通过 3N 单次传播获取 2N 数据点
disp('通过 3N 物理传播快速获取散点数据...');
num_runs = 5; % 运行 5 次独立的 3N 传播，即可产生 5 * 201 = 1005 个不重复的数据点！
H_data = [];
I_data = [];

for run = 1:num_runs
    % 随机产生独立的 S0 和 S_flip
    sigma = 2 * randi([0, 1], 1, n, 'gpuArray') - 1;   
    sigma_flip = 2 * randi([0, 1], 1, n, 'gpuArray') - 1;   
    
    % 单次 3N 物理传播，直接获取 2N+1 个时间片的光强(包含初态)
    I_windows = get_intensities_3N_single_pass(sigma, sigma_flip, Ugpu, R, TPIM, deltaf, dt, nt, idx_win0);
    
    % 理论精确重构这 2N+1 个滑动窗口对应的 spin 状态，并计算H
    for sel_idx = 1:(2*n+1)
        if sel_idx == 1
            sigma_candidate = sigma;
        else
            j_sel = sel_idx - 1;
            if j_sel <= n
                sigma_candidate = [sigma_flip(1 : j_sel), sigma(j_sel + 1 : n)];
            else
                k_sel = j_sel - n;
                sigma_candidate = [sigma(1 : k_sel), sigma_flip(k_sel + 1 : n)];
            end
        end
        
        % 记录数学 Hamiltonian(带有对角线的) 和对应的物理光强I
        H_math = -sigma_candidate * Jtest * sigma_candidate';
        H_data = [H_data, gather(H_math)];
        I_data = [I_data, gather(I_windows(sel_idx))];
    end
end

%% 4.  Fig.3 HI曲线散点图
figure(3);
fs = 22;    % 字体大小
scatter(H_data + Hdiag, I_data / I2E, 35, 'ob', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
% 绘制理想状态的参考虚线（用Hdiag来表示对角线的贡献）
plot([min(H_data) max(H_data)]+ Hdiag, [min(H_data) max(H_data)] , 'r--', 'LineWidth', 2);
xlabel('H/J', 'FontName', 'Times New Rome', 'FontSize', fs);
ylabel('-I/J', 'FontName', 'Times New Rome', 'FontSize', fs);
ylim([-700 0])
% xlim([-400 100])
set(gca, 'FontSize', fs); 
grid off;
box on;
hold off;
%% 存储数据
% save('HIJ11.mat',"J11","H_data",'I_data');
%% 3N 物理传播函数
function I_windows = get_intensities_3N_single_pass(sigma, sigma_flip, Ugpu, R, TPIM, deltaf, dt, nt, idx_win0)
    n = size(Ugpu, 1);
    % 3N 滑动窗口会产生 2N+1 个干涉光强采样点
    I_windows = zeros(1, 2*n+1, 'gpuArray');
    I_tot_waveform = zeros(1, nt, 'gpuArray');
    
    % 在时域拼接成 SSflipS (3N 长度) 的物理自旋状态列车
    total_sequence = [sigma, sigma_flip, sigma]; 
    
    % 产生平坦光梳谱
    freqa = ones(1, n, 'gpuArray');
    TPIM.create_combs(n, freqa);
    
    % 一次性计算 R 个解耦通道的 3N 物理调制与传播
    for rJ = 1:R
        xi_k = Ugpu(:, rJ)';                  
        xi_k_replicated = [xi_k, xi_k, xi_k]; % 同样进行 3N 长度复制
        
        % 调制信号
        initheav = total_sequence .* xi_k_replicated; 
        
        % 调制并传播
        % 调制时需要在时域第 K 块脉冲上加载相应的交替相位补偿信号 (-1)^K
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