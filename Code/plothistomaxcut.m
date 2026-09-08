%% maxcut histogram 绘制
clc;

% 1. 加载数据
data1 = Energy_minofmethod3N;
data2 = Energy_min_singlespin;

% 自动提取两组数据的能量值和频数
[counts1, centers1] = groupcounts(data1(:)); 
[counts2, centers2] = groupcounts(data2(:)); 

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

% 2. 开始画图
fig = figure('Color', 'w');
% 设置插图尺寸 (12cm x 9.5cm)
set(fig, 'Units', 'centimeters', 'Position', [5, 5, 12, 9.5]);
fs = 28;
% 颜色定义
color_tpim = [0.35, 0.42, 0.95]; % 蓝
color_sa   = [0.95, 0.45, 0.35]; % 橙红 (或 0.98, 0.95, 0.35 亮黄)

% 绘制柱状图
b1 = bar(unique_centers, y1, 'BarWidth', 1, ...
    'FaceColor', color_tpim, 'FaceAlpha', 0.65, 'EdgeColor', 'k', 'LineWidth', 0.8);
hold on;
b2 = bar(unique_centers, y2, 'BarWidth', 1, ...
    'FaceColor', color_sa, 'FaceAlpha', 0.65, 'EdgeColor', 'k', 'LineWidth', 0.8);

% X 轴刻度需要范围
x_min_val = min(unique_centers);
x_max_val = max(unique_centers);
span = x_max_val - x_min_val;

% 自动生成若干刻度，横向平铺
xticks(linspace(floor(x_min_val/10)*10, ceil(x_max_val/10)*10, 3));
set(gca, 'XTickLabelRotation', 0); 

% 设置坐标轴， X 轴左右两端留出适当边距
xlim([x_min_val - span*0.1, x_max_val + span*0.4]);
ylim([0, 100]);
box on;
grid off;
xlabel('Ising energy', 'FontName', 'Times New Roman', 'FontSize', fs);
ylabel('number of runs', 'FontName', 'Times New Roman', 'FontSize', fs);
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'LineWidth', 1.0);

% 图例 
lgd = legend([b1, b2], {'TPIM (algo)', 'SA'}, 'Location', 'northeast');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'Color', 'w', 'EdgeColor', 'k');

hold off;