%% 多组退火能量曲线绘制 

% 数据参数设置，iterations设置成实际的，例如***3040.mat表示30*40=1200
num_runs = 100;
num_iterations = 2000;
Hdiag = trace(Jtest);

% 1.1 可调参数
% 选择要绘制的实验组索引
runs_to_plot = [5,17,29,41,53,65,77,89]; 
M = length(runs_to_plot); % 绘制的退火组数

line_width = 1;   % 曲线的线宽
font_size = 28;   % 坐标轴和图例字体大小

% 1.2 配色
COLOR_PALETTE = [
    0.12 0.47 0.71;  % Blue (深蓝)
    0.83 0.37 0.00;  % Orange (橙红)
    0.30 0.69 0.29;  % Green (绿)
    0.78 0.10 0.78;  % Purple (紫)
    0.50 0.50 0.50;  % Gray (灰)
    0.90 0.60 0.00;  % Yellow-Orange (黄橙)
    0.20 0.80 0.80;  % Cyan (青)
    0.95 0.50 0.60;  % Pink (粉)
    0.60 0.40 0.20;  % Brown (棕)
    0.00 0.00 0.00   % Black (黑)
];

% 确保颜色数量足够
if M > size(COLOR_PALETTE, 1)
    warning('绘制的组数超过了预设的颜色数量，颜色可能会重复。');
    COLOR_PALETTE = lines(M); 
end


% 2: 绘制曲线
fig = figure('Color', 'w', 'Position', [150, 150, 630, 420]);
hold on;

legend_entries = cell(M, 1); % 存储图例标签

for idx = 1:M
    run_index = runs_to_plot(idx); % 获取要绘制的原始行索引
    
    % 获取能量数据(实际考虑的是分解模型，去掉对角线后代表实际哈密顿量)
    energy_data = Energy_iteration(run_index, :) + Hdiag;
    
    % 获取颜色
    current_color = COLOR_PALETTE(idx, :); 
    
    % 绘制曲线 
    % 绘制平滑曲线
    plot(1:num_iterations, energy_data/1000, ...
         'Color', current_color, ...
         'LineWidth', line_width);

    % 设置图例标签
    legend_entries{idx} = ['Run ', num2str(run_index)];
end

xlabel('Iterations', 'FontName', 'Times New Roman', 'FontSize', font_size);
ylabel('H (10^3/J)', 'FontName', 'Times New Roman', 'FontSize', font_size);
grid off;
box on;
set(gca, 'FontSize', font_size); 
hold off
ylim([-1.5 0.1]);
yyaxis right
plot(1:num_iterations, Tall, 'k');
ax = gca;
ax.YColor = 'k';
ylabel('T/J','Color','k');
ylim([0 80]);

% 图例
lgd = legend('Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5', 'Run 6', 'Run 7', 'Run 8', ...
    'NumColumns', 2, ...
    'FontSize', 21, ...
    'Location', 'northeast', ...
    'FontName', 'Times New Roman');
lgd.ItemTokenSize = [45, 10];
lgd.Position = [0.35, 0.58, 0.45, 0.3];
% 确保 X 轴和 Y 轴从数据最小值开始，可以根据实际情况调整
xlim([0, 2100]); 
hold off;