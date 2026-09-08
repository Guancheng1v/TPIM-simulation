%% 双 Y 轴算法演化对比图
clc; 
%% ==================== 1. 用户控制参数 ====================
run_idx   = 1;     % 默认提取第几组独立实验的数据进行可视化 (1 到 100)
%% 2. 从矩阵中提取第 run_idx 组的演化向量
if exist('Energy_singlespin', 'var') && exist('Energy_method3N', 'var')
    H_single_traj = Energy_singlespin(run_idx, :); 
    H_3N_traj     = Energy_method3N(run_idx, :);  
else
    error('请先在 MATLAB 工作区中加载包含 Energy_singlespin 和 Energy_method3N 的数据！');
end
totstep = length(H_3N_traj);
%% 3. 自动提取或计算 Wtotal 与基态能量 H_GS
if exist('Energy_minofmethod3N', 'var')
    H_GS = min(Energy_minofmethod3N);
else
    H_GS = min([H_single_traj, H_3N_traj]);
end
W_tot = W3_total;
%% 4. 通用精确换算公式：从 H 换算为 Cut Size C
% 公式: C = (2 * Wtotal - H) / 4
C_single_traj = (2 * W_tot - H_single_traj) / 4;
C_3N_traj     = (2 * W_tot - H_3N_traj) / 4;
C_GS          = (2 * W_tot - H_GS) / 4;

%% 5. 精密计算左右 Y 轴范围 (确保递增向量传给 ylim)
% H_bottom 较小 (更负), H_top 较大
H_bottom = H_GS - 20; 
H_max_val = max([H_single_traj, H_3N_traj]);
H_top = H_max_val + 20;

% C_top 数值较小, C_bottom 数值较大
C_top    = (2 * W_tot - H_top) / 4;
C_bottom = (2 * W_tot - H_bottom) / 4;

%% 6. 双 Y 轴对比图
fig = figure('Color', 'w', 'Position', [150, 150, 630, 420]);
fs = 28;
% =========================================================
% 左轴为 Ising Energy H, 右轴为 Cut Size C
yyaxis left
h1 = plot(H_single_traj/1e3, 'LineWidth', 1.2, 'Color', [0.85, 0.35, 0.25], 'LineStyle', '--'); hold on;
h2 = plot(H_3N_traj/1e3,     'LineWidth', 1.5, 'Color', [0.12, 0.53, 0.90], 'LineStyle', '-');

ylabel('Ising Energy H/10^3', 'FontName', 'Times New Roman', 'FontSize', fs);
set(gca, 'YColor', [0.2, 0.2, 0.2]); % 右轴刻度深灰色

% H 自下而上变大 (递增向量 [H_bottom, H_top])
ylim([H_bottom/1e3, H_top/1e3]); 

% 绘制 Ground State 红色虚线与文字
parCGSstr = ['(C=', num2str(C_GS), ')'];
% 用 text 函数放置文字位置
text(0.74 * totstep, H_GS/1e3 + 0.28, 'Max-Cut ', ... 
'Color', 'r', 'FontName', 'Times New Roman', 'FontSize', 18);
text(0.74 * totstep, H_GS/1e3 + 0.14, parCGSstr, ... 
'Color', 'r', 'FontName', 'Times New Roman', 'FontSize', 18);
yyaxis right
% 绘制红色虚线 
yline(C_GS, 'r--', 'LineWidth', 1.6);
ylabel('Max-Cut Size C', 'FontName', 'Times New Roman', 'FontSize', fs);
set(gca, 'YColor', [0.2, 0.2, 0.2]); % 右轴刻度深灰色
ylim([C_top, C_bottom]); 
set(gca, 'YDir', 'reverse'); % 翻转右轴


xlabel('Iterations', 'FontName', 'Times New Roman', 'FontSize', fs);
legend([h2, h1], {'TPIM(algo)', 'SA'}, ...
    'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 21);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs);
box on;
hold off;