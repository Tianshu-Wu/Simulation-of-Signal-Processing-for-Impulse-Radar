% 基于CFAR（恒定虚警率）的雷达目标检测算法
function [coordinate_tar, TH] = CFAR_DetectByAm(image, threshold, cfar_window, Pose)

% 输入参数：
%   image      : 输入的雷达幅度图像矩阵，尺寸为 [距离单元数, 方位单元数]
%   threshold  : 虚警率阈值（单位为dB），用于控制检测灵敏度
%   cfar_window: 结构体，定义CFAR窗口参数，包含以下字段：
%                - pazi_size : 方位维保护单元尺寸（单侧）
%                - prange_size: 距离维保护单元尺寸（单侧）
%                - mazi_size : 方位维检测单元尺寸（单侧）
%                - mrange_size: 距离维检测单元尺寸（单侧）
%   Pose       : 算法选择标志，0表示CA-CFAR，1表示OS-CFAR
% 输出参数：
%   coordinate_tar: 检测到的目标信息矩阵，每行格式为 [距离坐标, 方位坐标, 目标功率, 检测阈值]
%   TH          : 每个像素点的动态检测阈值矩阵，与输入图像同尺寸

%% 参数初始化
if nargin == 3          % 若输入参数为3个，默认使用CA-CFAR
    Pose = 0;           % Pose=0 表示CA-CFAR算法
end

TH = zeros(size(image));% 初始化阈值矩阵
[M, N] = size(image);   % 获取雷达图像的尺寸（M:距离单元数，N:方位单元数）
target_num = 0;         % 目标计数器
coordinate_tar = [];    % 目标信息存储矩阵

%% 窗口参数解析
% 保护窗口尺寸（排除目标邻近单元）
pazi_size = cfar_window.pazi_size;   % 方位保护单元尺寸
pran_size = cfar_window.prange_size; % 距离保护单元尺寸
% 检测窗口尺寸（用于计算背景噪声）
tazi_size = cfar_window.mazi_size;   % 方位检测单元尺寸
tran_size = cfar_window.mrange_size; % 距离检测单元尺寸

%% 计算窗口半径（转换为单侧宽度）
pazi_w = ceil((pazi_size-1)/2);  % 方位保护单元单侧宽度
pran_w = ceil((pran_size-1)/2);  % 距离保护单元单侧宽度
tazi_w = ceil((tazi_size-1)/2);  % 方位检测单元单侧宽度
tran_w = ceil((tran_size-1)/2);  % 距离检测单元单侧宽度

%% 生成窗口索引
% 保护窗口索引范围（相对于中心点）
pazi_indice = (-pazi_w:pazi_w);  % 方位维保护单元索引
pran_indice = (-pran_w:pran_w);  % 距离维保护单元索引
% 检测窗口索引范围（相对于中心点）
tazi_indice = (-tazi_w:tazi_w);  % 方位维检测单元索引
tran_indice = (-tran_w:tran_w);  % 距离维检测单元索引

%% 计算窗口中心位置（用于排除保护区域）
centre_a = (length(tazi_indice)-1)/2 + 1; % 方位维检测窗口中心
centre_r = (length(tran_indice)-1)/2 + 1; % 距离维检测窗口中心

%% 滑动窗口遍历图像
% 注意：边界区域不处理（避免越界）
for nn = tazi_w+1 : N-tazi_w      % 遍历方位单元（跳过边界）
    for mm = tran_w+1 : M-tran_w  % 遍历距离单元（跳过边界）
        %% 提取检测窗口区域
        detect_area = image(mm+tran_indice, nn+tazi_indice); % 当前检测窗口
        
        %% 屏蔽保护单元（置零）
        detect_area(centre_r+pran_indice, centre_a+pazi_indice) = 0; % 将保护区域置零
        
        %% 提取有效样本
        Sample = reshape(detect_area, [], 1);  % 将窗口展平为列向量
        Sample(Sample == 0) = [];              % 移除保护单元（零值）
        
        %% CFAR算法核心
        switch Pose
            case 0  % CA-CFAR: 计算背景单元的平均功率
                power_mean = mean(abs(Sample).^2);
            case 1  % OS-CFAR: 排序后选择特定百分比样本
                Sample_sorted = sort(Sample, 'ascend');
                k = 5 * round(length(Sample_sorted)/6); % 选择第5/6分位点
                power_mean = abs(Sample_sorted(k)).^2;   % 取对应位置的功率
        end
        
        %% 计算动态阈值
        THRESHOLD = 10^(threshold/10);      % 将dB转换为线性值
        TH(mm, nn) = power_mean * THRESHOLD; % 当前像素的检测阈值
        
        %% 目标判定
        Tcur = abs(image(mm, nn))^2;       % 当前像素的功率
        if Tcur >= TH(mm, nn)              % 超过阈值则判定为目标
            coordinate_tar = [coordinate_tar; mm, nn, Tcur, TH(mm, nn)];
            target_num = target_num + 1;   % 目标计数+1
        end
    end
end
end