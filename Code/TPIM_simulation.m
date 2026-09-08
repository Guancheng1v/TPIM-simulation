classdef TPIM_simulation < handle
    
    properties
        nt         % number of sample
        t_tol      % time limit
        time       % time linspace
        ff         % freq linspace
        
        lambda     % center wavelength
        deltaf     % interval of freq combs
        carrier    % combs in TD
        pulsenum   % pulse number in one rectangle modulation
        modu       % modulater in TD
        input      % signal before propagation
        
        D          % propagation dispersion
        Elt_new    % signal after propagation
        M
        y_mark     % select accurate time output signal 
    end
    
    methods(Static)
        function output_TD = dispersion_propagation(input_TD,lambda0,totL,D,num_of_iteration,~,ff)
            % 色散光纤中的光信号传输
            % 输入：
            %   input_TD            - 输入光信号
            %   lambda0             - 光信号中心波长(m)
            %   totL                - 传播距离(km)
            %   D                   - 色散参量(ps/nm/km)
            %   num_of_iteration    - 传播取点个数(可直接取1)
            %   nt                  - 时间窗口取点个数
            %   ff                  - 频域的linspace
            % 输出：
            %   output_TD           - 输出光信号

            c_const = 3e8;
            omega = 2*pi*c_const/lambda0; 
            Elf = fftshift(fft(fftshift(input_TD)));      % 输入信号的频域
            j = 1;
            L = totL/num_of_iteration  * j;           % 传播距离,km
            b0 = 2*pi/lambda0;
            b1 = 0;                                   % beta1,ps/km, use coordinate moving with signal
            b2 = - D * lambda0 / omega * 1e21;        % beta2,psps/km

            omg_i = 2*pi*ff;
            phi = (b0+b1*1e-12*(omg_i) + 0.5*b2*1e-24*(omg_i).^2)*L;
            Y_new = Elf.*exp(1i*phi);
            Elt_new = ifftshift(ifft(ifftshift(Y_new)));     % 时域的传播信号
            output_TD = Elt_new;
        end
    end
    
    methods
        function self = TPIM_simulation(nt0,t_tol0,deltaf0,lambda0,pulsenum0,D0)
            % 构造函数
            % 输入：
            %   nt0                 - 时间窗口取点个数
            %   t_tol0              - 时间窗口宽度
            %   deltaf0             - 频域梳齿间距(Hz)
            %   lambda0             - 光信号中心波长(m)
            %   pulsenum0           - 每个调制信号中的脉冲个数
            %   D0                  - 色散参量(ps/nm/km)
            %   num_of_iteration    - 传播取点个数(可直接取1)
            self.lambda = lambda0;
            self.t_tol = t_tol0;
            self.nt = nt0;
            dt = self.t_tol/(self.nt-1);
            df = 1/dt/self.nt;
            time_cpu = linspace(-self.t_tol/2, self.t_tol/2, self.nt);
            self.time = gpuArray(time_cpu);
            ff_cpu = (-self.nt/2+1:1:self.nt/2)*df;
            self.ff = gpuArray(ff_cpu);
            self.deltaf = deltaf0;
            self.pulsenum = pulsenum0;
            self.D = D0;      
        end
        
        function self = create_combs(self,freqnum0,freqa0)
            % 产生并绘制输入的脉冲梳齿
            % 输入：
            %   freqnum0            - 梳齿个数
            %   freqa               - 梳齿振幅大小
            df = self.ff(2) - self.ff(1);           % 频域相邻点间距
            carrierf = zeros(size(self.ff),'gpuArray');        % 频域信号
            deltafnum = round(self.deltaf/df);      % 相邻梳齿间隔对应的取点个数
            freqnum = freqnum0;                     
            freqa = gpuArray(freqa0);                         

            igpu = gpuArray(-floor(freqnum/2):1:(ceil(freqnum/2)-1));
            fdgpu = self.nt/2 + igpu*deltafnum + 1;
            carrierf(fdgpu) = freqa(igpu + floor(freqnum/2) + 1);
            self.carrier = ifftshift(ifft(ifftshift(carrierf)));    % time comb

%             figure(1)
%             set(gca,'FontName','Times New Rome','FontSize',20);
%             plot(self.ff/1e9,abs(carrierf));
%             xlabel('$frequency(\rm{GHz})$','Interpreter','latex')
%             ylabel('$Amplitude$','Interpreter','latex')
%             xlim([-250 250]);
        end       
        
        function self = create_modu_ampandphase(self,cou,initheav,devianum,~)
            % 产生并绘制调制信号(振幅相位编码)
            % 输入：
            %   cou            - 调制信号的个数
            %   initheav       - 调制信号的实部
            %   devianum       - 调制信号的平移数量，一般取0
            self.modu = zeros(1,self.nt,'gpuArray');
            deltat = 1/self.deltaf;
            % -- encode rectangles of input signal--
            signofheav = (-1).^(0:1:(cou-1));
            heav = initheav.*signofheav;

            i_vec = gpuArray(ceil(-cou/2):ceil(cou/2-1));
            j_vec = i_vec + devianum; 
            k_vec = i_vec - ceil(-cou/2) + 1; % heav 索引

            % 1. 计算所有脉冲的起始时间 T1 和结束时间 T2
            % Const 是固定的时间偏移量
            Const = -deltat * (ceil(self.pulsenum/2) - 0.5);
            
            T1_vec = Const + self.pulsenum * deltat * j_vec;    % 1 x cou
            T2_vec = Const + self.pulsenum * deltat * (j_vec + 1); % 1 x cou
            
            % 2. 构造时间和阈值矩阵（cou x nt）
            Time_Matrix = repmat(self.time, cou, 1);    
            T1_Matrix = repmat(T1_vec, self.nt, 1)'; 
            T2_Matrix = repmat(T2_vec, self.nt, 1)'; 

            % 3. 向量化阶跃函数
            % Pulse_Matrix (cou x nt): 每一行是一个矩形脉冲 (1或0)
            Pulse_Matrix = double(Time_Matrix >= T1_Matrix) - double(Time_Matrix >= T2_Matrix);

            % 4. 获取幅度向量并构造幅度矩阵
            Heav_Vector = heav(k_vec)'; % cou x 1
            Heav_Matrix = repmat(Heav_Vector, 1, self.nt); % cou x nt

            % 5. 权重叠加求和
            % self.modu (1 x nt) 是所有脉冲的加权和
            Weighted_Pulse_Matrix = Heav_Matrix .* Pulse_Matrix; % cou x nt
            self.modu = sum(Weighted_Pulse_Matrix, 1);
               
%             figure(2);
%             plot(self.time/1e-12, real(self.modu))
%             set(gca,'FontName','Times New Rome','FontSize',12);
%             xlabel('$time (\rm{ps})$','Interpreter','latex')
%             ylabel('$Real (a.u.)$',Interpreter='latex')
%             xlim([-200 200]);
%             title('amplitude and phase modulate-TD');
%             grid on                                           
        end
        
        function self = create_inputsignal(self)
            % 产生传播前信号
            self.input = self.carrier.*self.modu;
%             Elf = fftshift(fft(fftshift(self.input)));
%             figure(3);
%             subplot(2,1,1);
%             plot(self.time/1e-12, abs(self.input))
%             title('After modulate-TD');
%             set(gca,'FontName','Times New Rome','FontSize',12);
%             xlabel('$time (\rm{ps})$','Interpreter','latex')
%             ylabel('$E (a.u.)$',Interpreter='latex')
%             xlim([-200 200]);
%             ylim([0 1e-3])
%             grid on
% 
%             subplot(2,1,2);
%             plot(self.ff/1e12, abs(Elf), 'LineWidth',2);
%             title('After modulate-FD')
%             set(gca,'FontName','Times New Rome','FontSize',12);
%             xlabel('$\omega/2\pi (\mathrm{THz})$',Interpreter='latex');
%             ylabel('$E(\omega)(a.u.)$',Interpreter='latex');
%             grid on
%             xlim([-0.5 0.5]);
        end
        
        function self = create_propagation(self,D,num_of_iteration)
            % 产生传播后的时域信号
            % 输入：
            %   D                      - 色散系数(ps/nm/km)
            %   num_of_iteration       - 传播迭代的个数
            c_const = 3e8;
            Deltat = self.pulsenum / self.deltaf;                       % 矩形调制窗口的大小(s)
            totL = Deltat * 1e12 / self.D /(self.deltaf * self.lambda^2 / c_const * 1e9);                  % 传播距离
            self.Elt_new = self.dispersion_propagation(self.input,self.lambda,totL,D,num_of_iteration,self.nt,self.ff);         
%             figure(4)
%             plot(gather(self.time/1e-12), gather(abs(self.Elt_new)))
%             set(gca,'FontName','Times New Rome','FontSize',12);
%             xlabel('$time (\rm{ps})$','Interpreter','latex')
%             ylabel('$E (a.u.)$',Interpreter='latex')
%             xlim([-400 400]);
%             title(['After prop' , num2str(totL) , 'km-TD']);
        end
        
        function self = plot_output_withphase(self)
            % 绘制传播后的时域信号(带相位)和输入的调制信号(振幅)
            
            % 获取当前图形窗口的大小
            fig = figure(7);
            current_pos = get(fig, 'Position');
            % 调整宽度（拉长图形）
            new_width = current_pos(3) * 1.7; 
            set(fig, 'Position', [current_pos(1), current_pos(2), new_width, current_pos(4)]);
            set(gca,'FontName','Times New Rome','FontSize',18);
            yyaxis left
            % 绘制归一化的信号振幅和相位
            plot(gather(self.time/1e-12), gather(abs(self.Elt_new)/max(abs(self.Elt_new))),'r')
            xlim([-75 75]);
            xlabel('t(ps)')
            ylim([-2 2]);
            phase = angle(self.Elt_new);   % 相位，范围 [-π, π]
            % 将相位归一化到 [0, 1] 用于 HSV 颜色映射
            hue = (phase + pi) / (2 * pi); % 归一化相位
            % 创建 HSV 颜色映射
            cmap = hsv(256);                                            % 生成 HSV 颜色映射
            color_indices = round(hue * (size(cmap, 1) - 1)) + 1;       % 将相位映射到颜色索引
            colors = cmap(color_indices, :);                            % 获取对应的 HSV 颜色
            
            hold on;
            for i = 1:(length(self.time)-1)
                % 填充颜色区域
                fill([self.time(i)/1e-12, self.time(i+1)/1e-12, self.time(i+1)/1e-12, self.time(i)/1e-12], ...
                     [0, 0, abs(self.Elt_new(i))./max(abs(self.Elt_new)), abs(self.Elt_new(i + 1))./max(abs(self.Elt_new))], ...
                     colors(i, :), 'EdgeColor', 'none');
            end
            colormap(hsv);                                                      % 设置颜色映射为 HSV
            colorbar('Ticks', [0, 0.5, 1], 'TickLabels', {'-\pi', '0', '\pi'}); % 添加颜色条
            hold off
            set(gca,'FontName','Times New Rome','FontSize',18);
            ylabel('Output Amplitude')

            yyaxis right
            % 绘制输入的矩形调制
            plot(self.time/1e-12, real(self.modu)./max(abs(self.modu)),'--','Color',[239/256,73/256,104/256],'LineWidth',2)
            ylabel('Input \sigma')
            ylim([-1.1 1.1])
            set(gca, 'XAxisLocation', 'origin');
        end  
    end
    
end
