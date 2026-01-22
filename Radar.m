clc;clear;close all


%% 雷达与目标参数设置
c = 3e8; % 光速(m/s)

fc = 1e9;      % 载波频率(Hz)
lambda = c/fc; % 载波波长(m)
PRT = 800e-6; % 脉冲重复间隔(秒)
PRF = 1/PRT;  % 脉冲重复频率(Hz)
TW = 160e-6;  % 脉冲宽度(秒)
BW = 1e6;     % 带宽(Hz)
gamma = BW/TW; % 调频斜率
Ntm = 100;     % 脉冲积累个数

R_table = [25e3 20e3 20e3]; % 三个目标的距离(m)
Vs_table = [50 0 30];       % 三个目标的速度(m/s)



%% 生成线性调频信号回波
sif = 0; % 初始化中频信号
% 三个目标的回波信号叠加
for iT = 1:length(R_table)
    At = 1; % 目标回波幅度
    tm = (-Ntm/2:Ntm/2-1)*PRT; % 慢时间轴(脉冲间时间)
    R = R_table(iT);  % 当前目标初始距离
    Vs = Vs_table(iT); % 当前目标速度
    Rtm = R + Vs*tm;  % 每个脉冲对应的目标距离(考虑运动)
    
    Tt = 2*Rtm/c; % 目标回波延时(双程时间)
    
    frf = c/lambda; % 载波频率(Hz)
    fif = 60e6;     % 中频频率(Hz)
    fs = 4/3*fif;   % 中频采样率(Hz)
    
    R1 = 5e3;  % 近距门距离(m)
    R2 = 60e3; % 远距门距离(m)
    NPRT = 2*(R2-R1)/c * fs; % 快时间采样点数
    that = (0:NPRT-1)/fs + R1*2/c; % 快时间采样轴(脉冲内时间)
    
    [THAT TT] = meshgrid(that,Tt); % 创建网格矩阵用于计算
    
    % 中频回波信号模型（下变频后）:
    % 1. rectpulsef: 矩形窗函数，表示脉冲存在的时间范围
    % 2. cos项: 包含中频、载频和多普勒相位
    % 3. 最后一项: 线性调频信号的二次相位
    sif = sif + At .* rectpulsef(THAT-TT,TW) .* cos(2*pi.*fif.*THAT -2*pi*frf.*TT + pi*gamma.* (THAT-TT).^2);
end

% 绘制回波时域波形
figure;
imagesc(that*c/2*1e-3,tm,sif); % 距离转换为km，时间轴为慢时间
xlabel('距离 (km)');ylabel('慢时间 (s)');title('回波时域波形');
colorbar;



%% 加入噪声
Ps = At^2/2; % 信号功率
SNR = -10;    % 信噪比(dB)
Pn = Ps * 10^(-SNR/10); % 噪声功率
rng(100);    % 设置随机数种子保证结果可重复
noise = sqrt(Pn) * randn(size(sif)); % 生成高斯白噪声
sif = sif + noise; % 添加噪声后的中频信号

% 绘制第一个脉冲的波形
figure;
plot(that*c/2*1e-3,sif(1,:)); 
xlabel('距离 (km)');
ylabel('幅度');
title('第一个脉冲的波形(含噪声)');



%% 混频-中频正交采样
% 正交解调: 将中频信号下变频到基带
x_I = sif .* cos(2*pi*fif*THAT); % I路(同相分量)
x_Q = -sif .* sin(2*pi*fif*THAT); % Q路(正交分量)

% 绘制Q路信号波形
figure;
plot(x_Q(1,:));title('Q路信号波形');xlabel('采样点');ylabel('幅度');



%% 低通滤波——降采样
%正交采样滤波器设计
% 使用Parks McClellan算法设计FIR低通滤波器
fl = 1000; % 滤波器阶数
% 滤波器频带边缘: [0 过渡带起点 过渡带终点 1]
fbe = [0 (BW/2)/(fs/2) (BW/2+BW/5)/(fs/2) 1]; 
Damps = [1 1 0 0]; % 期望幅频响应
b = firpm(fl, fbe, Damps); % 设计滤波器系数

% 绘制滤波器频率响应
figure;
freqz(b);title('FIR滤波器的幅频和相频特性');

%滤波处理
% 对I/Q两路信号进行滤波
x_I_filter = filter(b, 1, x_I, [], 2); % 沿快时间维度(第2维)滤波
x_I_filter = x_I_filter(:, fl/2+1:end); % 去除滤波器延迟

x_Q_filter = filter(b, 1, x_Q, [], 2); % Q路滤波
x_Q_filter = x_Q_filter(:, fl/2+1:end); % 去除延迟

% 绘制滤波后的I/Q信号
figure;
subplot(2,1,1); 
plot(x_I_filter(1,:));ylabel('x_I_filter');title('低通滤波后的I路信号');

subplot(2,1,2); 
plot(x_Q_filter(1,:));ylabel('x_Q_filter');title('低通滤波后的Q路信号');

% 绘制Q路信号频谱
Nfft = size(x_Q_filter, 2);
freq = (-Nfft/2:Nfft/2-1)/Nfft*fs;
figure;
plot(freq, abs(fftshift(fft(x_Q_filter(1,:)))));title('Q路信号频谱');
xlabel('频率 (Hz)');ylabel('幅度');xlim([-.5e7,.5e7]);

%抽样(降采样)
x_I_ex = x_I_filter(:, 1:2:end); % I路降采样(每隔一个点取一个)
x_Q_ex = x_Q_filter(:, 1:2:end); % Q路降采样

% 绘制抽样后的I/Q信号
figure;
subplot(2,1,1);
plot(x_I_ex(1,:));ylabel('x_I_ex');title('抽样后的I路信号');

subplot(2,1,2);
plot(x_Q_ex(1,:));ylabel('x_Q_ex');title('抽样后的Q路信号');

% 合成复信号
x_ex = x_I_ex + 1i*x_Q_ex; % I路为实部，Q路为虚部

% 绘制合成信号
figure;
plot(real(x_ex(1,:)));ylabel('x_ex');title('正交下变频的最终合成信号(实部)');

% 生成降采样后的线性调频信号
Np = TW * fs/2; % 波形宽度对应的点数
that1 = (-Np/2:Np/2-1)/(fs/2); % 时间轴
sig_t = exp(1i*pi*gamma*that1.^2);

% 绘制LFM信号实部
figure;
plot(that1, real(sig_t));
title('LFM信号实部波形');xlabel('时间 (s)');ylabel('幅度');

% 绘制LFM信号频谱
spec_sig = fftshift(fft(sig_t));
Nfft = length(that1);
freq = (-Nfft/2:Nfft/2-1)/Nfft*fs/2;
figure;
plot(freq, abs(spec_sig));
title('LFM信号频谱');xlabel('频率 (Hz)');ylabel('幅度');
xlim([-0.2e7,0.2e7]);



%% 脉冲压缩
% 脉冲压缩处理
h = conj(fliplr(exp(1i*pi*gamma*that1.^2))); % 匹配滤波器系数(时间反转的发射信号共轭)
Nfft = length(h) + size(x_ex, 2) - 1; % FFT长度

% 频域脉冲压缩(匹配滤波)
spc = ifft((ones(Ntm,1) * fft(h, Nfft)) .* fft(x_ex, Nfft, 2), Nfft, 2);
spc = spc(:, length(h)/2+1:end); % 去除暂态点

% 绘制脉冲压缩结果
NPRT = 2*(R2-R1)/c * fs/2; % 快时间采样点数
that = (0:NPRT-1)/(fs/2) + R1*2/c; % 时间轴
spec2 = spc(:, 1:fix(NPRT)); % 取有效部分

figure;
imagesc(that*c/2*1e-3, tm, abs(spec2));
xlabel('距离 (km)');ylabel('慢时间 (s)');title('脉冲压缩结果');
colorbar;


%% MTI与MTD
% MTI(动目标显示)处理
spec3 = spec2(2:end,:) - spec2(1:end-1,:); % 单脉冲对消(消除静止目标)

% 绘制MTI结果
figure;
imagesc(that*c/2*1e-3, tm, abs(spec3));
title('MTI结果');xlabel('距离门');ylabel('脉冲数');
colorbar;

% MTD(动目标检测)处理
spec_coh = fftshift(fft(spec3, [], 1), 1); % 沿慢时间维做FFT(多普勒处理)
Ntm2 = Ntm-1; % 单脉冲对消后脉冲数减少1
fd = (-Ntm2/2:Ntm2/2-1)/Ntm2 * PRF; % 多普勒频率轴
vel = -fd*lambda/2; % 将多普勒频率转换为速度

% 绘制MTD结果(距离-速度三维图)
figure;
mesh(that*c/2*1e-3, vel, abs(spec_coh).^2);
title('MTD结果');xlabel('距离 (km)');ylabel('速度 (m/s)');zlabel('功率');



%% CA-CFAR(恒虚警率检测)处理
cfar_window.pazi_size = 10; % 方位向(多普勒维)保护单元数
cfar_window.prange_size = 10; % 距离向保护单元数
cfar_window.mazi_size = 20; % 方位向参考单元数
cfar_window.mrange_size = 20; % 距离向参考单元数
threshold = 15; % 检测门限(dB)

% 执行CFAR检测
[coordinate_tar, TH] = CFAR_DetectByAm(abs(spec_coh).^2, threshold, cfar_window, 0);
TH(find(TH==0)) = nan; % 将零值替换为NaN便于绘图

% 绘制MTD结果(含CFAR门限)
figure;
mesh(that*c/2*1e-3, vel, abs(spec_coh).^2);
hold on;
mesh(that*c/2*1e-3, vel, abs(TH));
title('MTD结果(含CFAR阈值门限)');
xlabel('距离 (km)');ylabel('速度 (m/s)');zlabel('功率');

% 绘制特定多普勒通道的距离剖面
figure;
plot(that*c/2*1e-3, abs(spec_coh(71,:)));
xlabel('距离 (km)');ylabel('幅度');title('MTD速度切片(第71通道)');