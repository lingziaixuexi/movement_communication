function h = jakes_rayleigh(N, fs, fd, M)
%JAKES_RAYLEIGH 改进Jakes(Zheng-Xiao)正弦波叠加法生成瑞利衰落信道
%   h = jakes_rayleigh(N, fs, fd, M)
%   N  : 样本数
%   fs : 采样率 (Hz)
%   fd : 最大多普勒频移 (Hz)
%   M  : 正弦波叠加支路数 (默认 32)
%   输出 h 为 N x 1 复信道增益, E[|h|^2] = 1
%
%   模型: Y.R. Zheng & C. Xiao, "Improved models for the generation of
%   multiple uncorrelated Rayleigh fading waveforms," IEEE Commun. Lett., 2002.

if nargin < 4
    M = 32;
end

t  = (0:N-1).' / fs;
hi = zeros(N, 1);
hq = zeros(N, 1);
theta = (2*rand - 1) * pi;          % 随机到达角偏置
for n = 1:M
    alpha = (2*pi*n - pi + theta) / (4*M);   % 第n条径的到达角
    phi1  = (2*rand - 1) * pi;               % 同相支路随机初相
    phi2  = (2*rand - 1) * pi;               % 正交支路随机初相
    hi = hi + cos(2*pi*fd*t*cos(alpha) + phi1);
    hq = hq + cos(2*pi*fd*t*sin(alpha) + phi2);
end
% 归一化: 每支路方差 M/2, 除以 sqrt(M) 后 I/Q 各 1/2, 总功率为 1
h = (hi + 1j*hq) / sqrt(M);
end
