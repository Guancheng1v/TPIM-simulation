%% cholesky 分解
% n为自旋的个数，J为相互作用矩阵 
function [Ksi,Jtest]=JtoKsi(n,J,c)
    TOL = 1e-9; % 引入数值容差
    if nargin == 3
        for i=1:1:n
            J(i,i) = c;
        end
    end
    Ksi=zeros(n,n);
    i=n;
    while(i>0)
        ksi=zeros(i,1);
        current_diag = J(1,1); % 当前子矩阵的对角线元素 J_ii
        % 1. 检查 J(1,1) 是否接近或小于零
        if current_diag < TOL 
            % 如果 J(1,1) 接近或小于零，则认为分解已经完成 (秩 K 已用完)
            ksi(1,1) = 0;
            ksi(2:end,1) = 0; % 确保这一列的其余元素也为零
        else
            % 2. 标准 Cholesky 分解步骤 (下三角)
            ksi(1,1) = sqrt(current_diag); 
            
            % 3. 计算 ksi 的其余元素 (除法步骤)
            ksi(2:end,1)=(J(1,2:end)/ksi(1,1)).';
        end
    
        % 4. 更新残差矩阵 J0
        J0=J-ksi*ksi.'; 
        
        % 5. 准备下一轮循环的子矩阵
        J=J0(2:end,2:end);
        
        % 6. 将当前计算的 ksi 向量放入最终的 Ksi 矩阵
        Ksi(n+1-i:end,n+1-i)=ksi;
        
        i=i-1;
    end 
    Jtest=Ksi*Ksi.';
end

