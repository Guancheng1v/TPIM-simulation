%% Ising histogram
fs = 28;
load('TPIM3N_multi_data_Jrank13040.mat')
% 可调右侧空白（如 0.8 表示空出 0.8 个柱子宽度，1.2 表示 1.2 个，以此类推）
right_padding = 0.55; 
Hdiag = trace(Jtest);

for i = 1:1:100
    data = Energy_min + Hdiag; 
end
% 1. 自动提取两组数据的能量值和频数
[counts1, centers1] = groupcounts(data(:)); 
[counts2, centers2] = groupcounts(Energy_min_singlespin(:) + Hdiag); 

% 2. 统一并对齐两组格点
unique_centers = unique([centers1(:); centers2(:)]);
unique_centers = sort(unique_centers, 'ascend'); % 升序

y1 = zeros(length(unique_centers), 1);
y2 = zeros(length(unique_centers), 1);
for k = 1:length(unique_centers)
    idx1 = find(centers1 == unique_centers(k));
    if ~isempty(idx1), y1(k) = counts1(idx1); end
    
    idx2 = find(centers2 == unique_centers(k));
    if ~isempty(idx2), y2(k) = counts2(idx2); end
end

L = length(unique_centers); % 真实的有效自旋能量点个数
x_coords = 1:L;             % 纯数值坐标轴

% 3. 开始画图
fig = figure('Color', 'w');
set(fig, 'Units', 'centimeters', 'Position', [5, 5, 16, 9.5]); 
color_tpim = [0.35, 0.42, 0.95]; % 学术蓝
color_sa   = [0.95, 0.45, 0.35]; % 学术橙红 (或 0.98, 0.95, 0.35 亮黄)

% 绘制图层 1：光子仿真（蓝色）
b1 = bar(x_coords, y1, 0.3, 'FaceColor', color_tpim, 'FaceAlpha', 0.65);
hold on;
grid on;

% 绘制图层 2：模拟退火（橙红色）
b2 = bar(x_coords, y2, 0.3, 'FaceColor', color_sa, 'FaceAlpha', 0.65);

% 4. 手动设置 X 轴显示边界
xlim([0.7, L + right_padding]); 

xlabel('H/J', 'FontName', 'Times New Roman', 'FontSize', fs);
ylabel('number of runs', 'FontName', 'Times New Roman', 'FontSize', fs);
ylim([0 100]);

% 5. 贴上 X 轴标签
set(gca, 'XTick', x_coords, 'XTickLabel', string(unique_centers),'Fontsize',fs);

% 6. 图例
lgd = legend([b1, b2], {'TPIM (simu)', 'SA'}, 'Location', 'northeast');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'Color', 'w', 'EdgeColor', 'k');
hold off;