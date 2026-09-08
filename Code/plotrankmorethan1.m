%% rank more than 1 绘制
%% 1.1 参数设置 
J = Jrank5;
rankofJJ = rank(J);
figWidth = 16;         
barHeight = 1.2;       
gapHeight = 0.5;          
bottomSpace = 1.6;      % 稍微增加底部空间，防止 Spin Index 太挤
topSpace = 0.3;         
leftSpace = 1.5;          
rightSpace = 2;        
totalHeight = bottomSpace + topSpace + (rankofJJ * barHeight) + ((rankofJJ-1) * gapHeight);
                
%% 1.2. 数据源与奇偶维度转换（将自旋规模 N 转换为奇数以对齐光梳中心）
J_num = size(J, 1);
if mod(J_num, 2) == 0
    J1gpu = zeros(J_num+1, J_num+1, 'gpuArray');
    J1gpu(1:end-1, 1:end-1) = gpuArray(J);
else
    J1gpu = gpuArray(J);
end
Hdiag = trace(J);
n = size(J1gpu, 1);

%% 1.3. Cholesky 矩阵分解，
% Cholesky分解，输出Xi矩阵和实际考虑的相互作用矩阵Jtest
[U,Jtest] = JtoKsi(n,J1gpu,1);
R = rank(U);
column_sums = sum(abs(U)); 

% 使用一个微小的阈值（如 1e-12）来判断哪些列是有值的
is_nonzero = (column_sums > 1e-12); 

% 只提取有值的列。这会将 U 矩阵的尺寸从 101x101 压缩为 101xR！
U = U(:, is_nonzero);
Ugpu = gpuArray(U);
disp(['Cholesky 分解完成，系统解耦维度（秩）为: ', num2str(R)]);
%% 2. 定义离散颜色映射 (红-白-蓝)
% 对应关系: 
% 第一行 -> 最小值 (-1) -> 蓝色
% 第二行 -> 中间值 (0)  -> 白色
% 第三行 -> 最大值 (1)  -> 红色
my_cmap = [
    0, 0, 1;   % 蓝色 (Blue) - 对应 -1
    1, 1, 1;   % 白色 (White) - 对应 0
    1, 0, 0    % 红色 (Red) - 对应 1
];

%% 3. 创建画布
fig = figure('Color', 'w');
set(fig, 'Units', 'centimeters', 'Position', [5, 5, figWidth, totalHeight]);
fs = 28;  
%% 4. 绘制每一个 Rank
for i = 1:rankofJJ
    currentY = bottomSpace + (rankofJJ - i) * (barHeight + gapHeight);
    axPos = [leftSpace/figWidth, currentY/totalHeight, ...
             (figWidth-leftSpace-rightSpace)/figWidth, barHeight/totalHeight];
    
    ax = axes('Position', axPos);
    
    % 设置范围为 [-1.5, 1.5]，这样 -1, 0, 1 会分别落在 [ -1.5, -0.5], [-0.5, 0.5], [0.5, 1.5] 三个区间内
    imagesc(U(:, i)', [-1.5, 1.5]); 
    colormap(ax, my_cmap); % 应用三色映射
    
    % 外观设置
    set(ax, 'Box', 'on', 'LineWidth', 0.8, 'TickDir', 'out');
    set(ax, 'YTick', [], 'XTick', []); 
    
    yl = ylabel(['\xi_', num2str(i)], 'FontSize', fs, 'Rotation', 0, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    yl.Units = 'normalized'; 
    yl.Position(1) = -0.03; 
    yl.Position(2) = 0.5;

    if i == rankofJJ
        set(ax, 'XTick', [1, 100]);
        set(ax, 'TickLength', [0, 0]); 
        xl = xlabel('Spin Index', 'FontSize', fs, 'FontName', 'Times New Roman');
        xl.Position(2) = 1.75; 
    end
end
set(gca, 'FontSize', fs); % 设置坐标轴刻度字体
%% 5. 添加全局 Colorbar
% 隐藏坐标轴放置 Colorbar
cbAx = axes('Position', [leftSpace/figWidth , bottomSpace/totalHeight, ...
            (figWidth-leftSpace-rightSpace)/figWidth, ...
            (rankofJJ*barHeight + (rankofJJ-1)*gapHeight)/totalHeight], 'Visible', 'off');
        
clim(cbAx, [-1.5, 1.5]); % 统一色标
colormap(cbAx, my_cmap);
cb = colorbar('eastoutside');

% 强制 Colorbar 刻度在 -1, 0, 1
cb.Ticks = [-1, 0, 1]; 

% 调整 colorbar 的位置
cbPos = cb.Position;
cbPos(1) = (figWidth - rightSpace + 0.4) / figWidth; 
cb.Position = cbPos;
set(gca, 'FontSize', fs);