function J = generate_mattis_matrix(N, K)
% GENERATE_RANDOM_BINARY_RANK_4_J_MATRIX
%
% 构造 N x N 的相互作用矩阵 J，保证其秩为 K ，且其分解矩阵 L 
% 仅包含 {-1, 0, 1} 元素。

% 输入参数:
% N: 矩阵的维度 (自旋数量)，例如 100。
% K: 矩阵的秩 (默认选取 K=1)。

    if nargin < 2
        K = 1;
    end
    
    % 确保随机性可复现
    rng('shuffle'); % 使用时间作为种子，每次运行结果不同
    
    % 1. 定义目标下三角矩阵 L
    L = zeros(N, N);
    
    % 2. 定义 K 个周期模式，以确保 K 列向量的线性独立性
    % 使用质数周期 [2, 3, 5, 7] 来保证模式间的差异性，也可人为规定如1*100
%     pattern_periods = ones(1,100); 
    pattern_periods = [2, 3 , 5 , 7,11]; 
    if length(pattern_periods) < K
        error('模式数量不足 K。');
    end

    % 3. 填充 L 的前 K 列，且只填充下三角部分 {-1, 0, 1}
    for j = 1:K % 循环 K 列
        
        period = pattern_periods(j);
        
        % 填充 L 的第 j 列 (从第 j 行到第 N 行)
        for i = j:N
            % 3a. 结构化：
            if mod(i-j+1, period) == 1
                % 3b. 随机选择值
                L(i, j) = randi([-1, 1]);
            end
        end
    end
   
    % 4. 构造对称矩阵 J: J = L * L'
    % J 自动具有秩 K 且是半正定的。
    J = L * L';
    J = (J + J') / 2; % 确保对称
    
    disp(['J 是一个稠密、秩 ', num2str(K), ' 的相互作用矩阵 (L 元素为 {-1, 0, 1})。']);
    
    % 验证秩是否合理
    calculated_rank = rank(J); 
    disp(['计算的秩: ', num2str(calculated_rank)]);

    % 5. 输出边表信息
    filename = 'interaction_edges.txt';
    fid = fopen(filename, 'w');
    numofv = 0;
    J = J - diag(diag(J));
    % 遍历每一个顶点 i
    for i = 1:size(J, 1)
        % 找到与顶点 i 相连的所有顶点 j
        % 注意此处重复输出了 (i,j) 和 (j,i)
        for j = 1:size(J, 2)
            edgew = J(i, j);
            
            % 只有当相互作用不为 0 时才输出
            if edgew ~= 0
                % 每一行输出一条边：源  顶点 目标顶点 权重
                fprintf(fid, '%d %d %.4f\n', i, j, edgew);
                numofv = numofv+1;
            end
        end
    end
    
    fclose(fid);
    disp(['边表数据已成功导出至 ', filename]);
end

