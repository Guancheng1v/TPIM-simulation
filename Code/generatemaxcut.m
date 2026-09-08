%% Generate_3_Benchmark_Instances
% 生成 3 种 100-Spin 图 (莫比乌斯图、0/-1 随机图、1/-1 图)
clc;
N = 100; % 节点规模

%% ==================== 实例 1：Möbius Ladder  ====================
Wmobius = zeros(N, N);
for i = 1:N
    left_node  = mod(i - 2, N) + 1;
    right_node = mod(i, N) + 1;
    opp_node   = mod(i - 1 + N/2, N) + 1;
    Wmobius(i, left_node)  = 1;
    Wmobius(i, right_node) = 1;
    Wmobius(i, opp_node)   = 1;
end
Jmobius = -Wmobius; % Ising 矩阵
W1_total = sum(sum(triu(Wmobius, 1))); % 上三角权重和 = 150

% save('J_mobius_100.mat', 'Jmobius', 'Wmobius', 'W1_total', 'N');

disp(['[实例 1 生成成功] Möbius Ladder 100: W_total = ', num2str(W1_total), ' | 矩阵秩 = ', num2str(rank(Jmobius))]);

%% ==================== 实例 2：0/-1 随机无权 Max-Cut 图 (概率 0.5) ====================
Wrand_tri = rand(N, N) < 0.5; 
Wrand_tri = triu(Wrand_tri, 1);
Wrand = Wrand_tri + Wrand_tri'; % 对称化无权图 W_ij ∈ {0, 1}
Jrand = -Wrand;
W2_total = sum(Wrand_tri(:)); % 总边数 ≈ 2475

% save('J_unweighted_random_100.mat', 'Jrand', 'Wrand', 'W2_total', 'N');

disp(['[实例 2 生成成功] 0/-1 随机无权图 : W_total = ', num2str(W2_total), ' | 矩阵秩 = ', num2str(rank(Jrand))]);

%% ==================== 实例 3：1/-1 随机受挫 Max-Cut 图 (Biased, 概率 0.5) ====================
% 随机选择 +1 (反铁磁 W_ij=1, J_ij=-1) 或 -1 (铁磁 W_ij=-1, J_ij=1)
rand_mask = (rand(N, N) < 0.5) * 2 - 1; % 随机 +1 或 -1
Wrandbia_tri = triu(rand_mask, 1);
Wrandbia = Wrandbia_tri + Wrandbia_tri'; % 包含 +1 和 -1 的受挫矩阵
Jrandbia = -Wrandbia;
W3_total = sum(Wrandbia_tri(:)); % 权重叠加和

% save('J_biased_random_100.mat', 'Jrandbia', 'Wrandbia', 'W3_total', 'N');

disp(['[实例 3 生成成功] 1/-1 随机受挫图 : W_total = ', num2str(W3_total), ' | 矩阵秩 = ', num2str(rank(Jrandbia))]);
disp('==================================================');
