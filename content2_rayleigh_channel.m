%% 内容2: 瑞利衰落信道仿真与统计特性分析
% 采用改进Jakes(正弦波叠加)法生成瑞利平坦衰落信道, 分析:
%   (1) 包络随时间变化(深衰落现象)  (2) 包络分布 vs 理论瑞利PDF
%   (3) 相位分布 vs 均匀分布        (4) 自相关函数 vs 贝塞尔函数理论值
%   (5) 多普勒功率谱 vs Jakes经典U形谱
clear; close all; rng(2026, 'twister');

fc = 2e9;                 % 载波频率 2 GHz
v  = 120 / 3.6;           % 移动速度 120 km/h -> m/s
c  = 3e8;
fd = v / c * fc;          % 最大多普勒频移
fs = 100e3;               % 采样率 100 kHz
N  = 2^19;                % 样本数 (约5.2 s)
M  = 32;                  % 叠加正弦波支路数

fprintf('载波 %.1f GHz, 速度 %.0f km/h, 最大多普勒频移 fd = %.1f Hz\n', ...
        fc/1e9, v*3.6, fd);

h   = jakes_rayleigh(N, fs, fd, M);
env = abs(h);             % 包络
phs = angle(h);           % 相位
fprintf('信道平均功率 E[|h|^2] = %.4f (理论1)\n', mean(env.^2));

%% 1. 包络随时间变化曲线 (dB)
t = (0:N-1).'/fs;
figure('Position', [100 100 820 420]);
idx = t <= 0.5;                       % 显示前0.5 s
plot(t(idx)*1e3, 20*log10(env(idx)), 'b');
hold on; yline(0, 'r--', 'LineWidth', 1.2);
grid on; ylim([-40 10]);
xlabel('时间 (ms)'); ylabel('包络 20lg|h(t)| (dB)');
title(sprintf('瑞利衰落信道包络 (f_d = %.0f Hz)', fd));
legend('信道包络', '均方根电平 0 dB');
exportgraphics(gcf, '../figures/fig2_1_瑞利包络时变曲线.png', 'Resolution', 150);

%% 2. 包络分布 vs 理论瑞利PDF
% E[|h|^2]=1 时每维方差 s2=0.5, 理论 f(r) = (r/s2)exp(-r^2/(2*s2)) = 2r*exp(-r^2)
figure('Position', [100 100 760 420]);
histogram(env, 80, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9], ...
          'EdgeColor', 'none');
hold on;
r  = linspace(0, 3.5, 400);
s2 = 0.5;
plot(r, r/s2 .* exp(-r.^2/(2*s2)), 'r-', 'LineWidth', 1.8);
grid on;
xlabel('包络 r'); ylabel('概率密度 f(r)');
legend('仿真直方图', '理论瑞利PDF');
title('瑞利衰落包络概率密度');
exportgraphics(gcf, '../figures/fig2_2_瑞利包络分布.png', 'Resolution', 150);

%% 3. 相位分布 vs 均匀分布
figure('Position', [100 100 760 420]);
histogram(phs, 60, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9], ...
          'EdgeColor', 'none');
hold on; yline(1/(2*pi), 'r-', 'LineWidth', 1.8);
grid on; ylim([0 0.25]);
xlabel('相位 \theta (rad)'); ylabel('概率密度 f(\theta)');
legend('仿真直方图', '理论均匀分布 1/2\pi');
title('瑞利衰落相位概率密度');
exportgraphics(gcf, '../figures/fig2_3_瑞利相位分布.png', 'Resolution', 150);

%% 4. 自相关函数 vs J0贝塞尔理论 + 多普勒功率谱
maxlag = round(3/fd*fs);                       % 取约3个 1/fd
[acf, lags] = xcorr(h, maxlag, 'normalized');
acf  = real(acf(lags >= 0));
tau  = (0:maxlag).'/fs;
acf_theory = besselj(0, 2*pi*fd*tau);          % 理论: R(tau)=J0(2*pi*fd*tau)

[psd, f] = pwelch(h, hann(8192), 4096, 16384, fs, 'centered');
% 理论Jakes谱: S(f) = 1/(pi*fd*sqrt(1-(f/fd)^2)), |f|<fd
fj = linspace(-0.999*fd, 0.999*fd, 600);
Sj = 1 ./ (pi*fd*sqrt(1 - (fj/fd).^2));

figure('Position', [100 100 980 400]);
subplot(1,2,1);
plot(tau*fd, acf, 'b', 'LineWidth', 1.2); hold on;
plot(tau*fd, acf_theory, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('归一化时延 f_d\tau'); ylabel('归一化自相关');
legend('仿真', '理论 J_0(2\pi f_d\tau)');
title('信道自相关函数');
subplot(1,2,2);
plot(f, 10*log10(psd), 'b'); hold on;
plot(fj, 10*log10(Sj), 'r--', 'LineWidth', 1.5);
grid on; xlim([-3*fd 3*fd]); ylim([-60 -10]);
xlabel('频率 (Hz)'); ylabel('PSD (dB/Hz)');
legend('Welch估计', '理论Jakes谱', 'Location', 'south');
title('多普勒功率谱');
exportgraphics(gcf, '../figures/fig2_4_自相关与多普勒谱.png', 'Resolution', 150);

%% 5. 不同多普勒频移下包络对比 (慢衰落 vs 快衰落)
fd_list = [10 100 500];
figure('Position', [100 100 820 560]);
for ii = 1:numel(fd_list)
    hk = jakes_rayleigh(round(0.5*fs), fs, fd_list(ii), M);
    tk = (0:numel(hk)-1).'/fs;
    subplot(3,1,ii);
    plot(tk*1e3, 20*log10(abs(hk)), 'b');
    grid on; ylim([-40 10]);
    ylabel('包络 (dB)');
    title(sprintf('f_d = %d Hz', fd_list(ii)));
end
xlabel('时间 (ms)');
exportgraphics(gcf, '../figures/fig2_5_不同多普勒包络对比.png', 'Resolution', 150);

fprintf('内容2完成, 图片已保存到 figures/ 目录\n');
