%% rank only 1 绘制
%% 1.1. 参数设置
J = Jrank1;
rankofJJ = rank(J);
displayRank = max(rankofJJ, 2); 

figWidth = 16;          
barHeight = 1.2;       
gapHeight = 1.6;          
bottomSpace = 1.6;      % 稍微增加底部空间，防止 Spin Index 太挤
topSpace = 0.3;         
leftSpace = 1.5;          
rightSpace = 2;       % 显著增加右侧留白，为"100"指示和 Colorbar 腾空间

barAreaHeightPhys = (displayRank * barHeight) + ((displayRank-1) * gapHeight);
totalHeight = bottomSpace + topSpace + barAreaHeightPhys;

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
%% 2. 定义颜色映射
my_cmap = [0, 0, 1; 1, 1, 1; 1, 0, 0];

%% 3. 创建画布
fig = figure('Color', 'w');
set(fig, 'Units', 'centimeters', 'Position', [5, 5, figWidth, totalHeight]);
fs = 28; 
%% 4. 绘制数据条
actualBarsHeightPhys = (rankofJJ * barHeight) + ((rankofJJ-1) * gapHeight);
verticalOffset = (barAreaHeightPhys - actualBarsHeightPhys) / 2;

for i = 1:rankofJJ
    currentY = bottomSpace + verticalOffset + (rankofJJ - i) * (barHeight + gapHeight);
    
    % axPos 的宽度会因为 rightSpace 变大而自动收缩，给右侧留出空间
    axPos = [leftSpace/figWidth, currentY/totalHeight, ...
             (figWidth-leftSpace-rightSpace)/figWidth, barHeight/totalHeight];
    
    ax = axes('Position', axPos);
    imagesc(U(:, i)', [-1.5, 1.5]);
    colormap(ax, my_cmap);
    
    set(ax, 'Box', 'on', 'LineWidth', 0.8, 'TickDir', 'out');
    set(ax, 'YTick', [], 'XTick', []); 
    
    yl = ylabel(['\xi_', num2str(i)], 'FontSize', fs, 'Rotation', 0, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    yl.Units = 'normalized'; 
    yl.Position = [-0.03, 0.5, 0];

    if i == rankofJJ
        % 调整 X 轴刻度标签的对齐方式
        set(ax, 'XTick', [1, 94]);
        set(ax, 'XTicklabel', [1, 100]);
        set(ax, 'TickLength', [0, 0]); 
        ax.XAxis.TickLabelInterpreter = 'tex';
        xl = xlabel('Spin Index', 'FontSize', fs, 'FontName', 'Times New Roman');
        % 调整 xlabel 的垂直高度，避免和数字重叠
        xl.Position(2) = 1.8; 
    end
end
set(gca, 'FontSize', fs);
%% 5. 添加全局 Colorbar
cbAx = axes('Position', [leftSpace/figWidth , bottomSpace/totalHeight, ...
            (figWidth-leftSpace-rightSpace)/figWidth, barAreaHeightPhys/totalHeight], 'Visible', 'off');
        
clim(cbAx, [-1.5, 1.5]);
colormap(cbAx, my_cmap);
cb = colorbar('eastoutside');

cb.Ticks = [-1, 0, 1]; 
cb.TickLabels = {'-1', '0', '1'}; 

% 将 Colorbar 放在轴右侧
cbPos = cb.Position;
cbPos(1) = (figWidth - rightSpace + 0.4) / figWidth; 
cbPos(3) = 0.5 / figWidth; % 保持 colorbar 宽度固定(约0.5cm)
cb.Position = cbPos;

set(gca, 'FontSize', fs);