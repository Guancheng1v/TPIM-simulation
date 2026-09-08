% === 绘制相互作用圆环图 ===

J = J33; % 读取当前工作区的主相互作用矩阵
N = size(J, 1);
num_visible = 32; % 截取展现的自旋数量，可选择更少或更多
half_vis = num_visible / 2;

% 1. 映射自旋的索引（截取最后 32 个自旋进行局部可视化）
visible_indices = (N - half_vis*2 + 1):N;

% 2. 提取子矩阵，并去除对角线自耦合
J_sub = J(visible_indices, visible_indices);
J_sub = J_sub - diag(diag(J_sub)); 

% 调整数据为 CPU 模式
J_sub = gather(J_sub);

% 3. 建立无向图（使用绝对值作为边的粗细权重）
G = graph(abs(J_sub), 'upper');

% 4. 计算圆环格点坐标（多出来的点表示省略号）
total_slots = 36; 
angles = linspace(0, 2*pi, total_slots + 1);
angles(end) = []; 

x_nodes = cos(angles(1:num_visible));
y_nodes = sin(angles(1:num_visible));
x_ellipsis = cos(angles(num_visible+1 : total_slots-1));
y_ellipsis = sin(angles(num_visible+1 : total_slots-1));

% 5. 提取每条边的原始正负号，并赋予对应的红/蓝色
num_edges = numedges(G);
edge_colors = zeros(num_edges, 3); % 初始化 RGB 颜色矩阵
edge_nodes = G.Edges.EndNodes;     % 获取每条边的两端节点序号

for e = 1:num_edges
    u = edge_nodes(e, 1);
    v = edge_nodes(e, 2);
    original_val = J_sub(u, v);    % 获取带符号的原始相互作用值
    
    if original_val > 0
        % 正耦合（铁磁相互作用）：学术红
        edge_colors(e, :) = [230/255, 75/255, 75/255]; 
    else
        % 负耦合（反铁磁相互作用）：学术蓝
        edge_colors(e, :) = [75/255, 110/255, 220/255]; 
    end
end

% 6. 如果线条有不同粗细，非线性设置线宽 (范围设为 2 ~ 5)；或者统一为相同粗细
if num_edges > 0
    weights = G.Edges.Weight;
    min_w = 2; % 最弱相互作用对应的线条粗细
    max_w = 5; % 最强相互作用对应的线条粗细
    if max(weights) == min(weights)
        edge_widths = ones(size(weights)) * min_w;
    else
        % 采用开方压缩，拉近强弱线条的视觉差距
        norm_weights = (weights - min(weights)) / (max(weights) - min(weights));
        edge_widths = min_w + sqrt(norm_weights) * (max_w - min_w);
    end
else
    edge_widths = 2;
end

% 7. 作图
figure('Color', 'w', 'Position', [200, 200, 500, 500]);
hold on;
p = plot(G, 'XData', x_nodes, 'YData', y_nodes, 'NodeLabel', {});
p.Marker = 'o';
p.MarkerSize = 5.5;                     % 稍微放大节点
p.NodeColor = [70/255, 70/255, 70/255]; % 炭灰色

if num_edges > 0
    p.LineWidth = edge_widths;       % 赋予线宽
    p.EdgeColor = edge_colors;       % 赋予颜色
    p.EdgeAlpha = 0.5;               
end

% 8. 绘制代表未展现自旋的省略号 "..." (与节点颜色保持一致)
scatter(x_ellipsis, y_ellipsis, 14, [70/255, 70/255, 70/255], 'filled');

xlim([-1.08 1.08]);
ylim([-1.08 1.08]);
axis equal;
axis off; 
hold off;