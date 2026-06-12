%% 内容1: 高斯噪声的产生与统计特性验证
% 产生均值为0、方差为sigma^2的实高斯白噪声以及复高斯白噪声,
% 从时域波形、幅度概率密度、自相关函数、功率谱密度四个角度验证其统计特性。
clear; close all; rng(2026, 'twister');

fs     = 1e5;          % 采样率 (Hz)
N      = 2^18;         % 样本数
sigma2 = 1;            % 目标方差
sigma  = sqrt(sigma2);

%% 1. 产生实高斯白噪声
n_real = sigma * randn(N, 1);

fprintf('实高斯噪声: 均值 = %.4f (理论0), 方差 = %.4f (理论%.1f)\n', ...
        mean(n_real), var(n_real), sigma2);

% 复基带高斯白噪声: 总方差 sigma2, I/Q 各 sigma2/2
n_cplx = sqrt(sigma2/2) * (randn(N,1) + 1j*randn(N,1));
fprintf('复高斯噪声: 均值 = %.4f%+.4fj, 总方差 = %.4f (理论%.1f)\n', ...
        real(mean(n_cplx)), imag(mean(n_cplx)), var(n_cplx), sigma2);

%% 2. 时域波形
t = (0:N-1).'/fs;
figure('Position', [100 100 760 420]);
plot(t(1:1000)*1e3, n_real(1:1000), 'b');
grid on;
xlabel('时间 (ms)'); ylabel('幅度');
title('高斯白噪声时域波形 (前1000个样本)');
exportgraphics(gcf, '../figures/fig1_1_噪声时域波形.png', 'Resolution', 150);

%% 3. 幅度分布直方图与理论PDF对比
figure('Position', [100 100 760 420]);
histogram(n_real, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9], ...
          'EdgeColor', 'none');
hold on;
x = linspace(-5*sigma, 5*sigma, 500);
pdf_theory = 1/(sqrt(2*pi)*sigma) * exp(-x.^2/(2*sigma2));
plot(x, pdf_theory, 'r-', 'LineWidth', 1.8);
grid on;
xlabel('幅度 x'); ylabel('概率密度 f(x)');
legend('仿真直方图', '理论高斯PDF', 'Location', 'northeast');
title(sprintf('高斯噪声幅度分布 (\\mu=0, \\sigma^2=%.1f)', sigma2));
exportgraphics(gcf, '../figures/fig1_2_噪声幅度分布.png', 'Resolution', 150);

%% 4. 自相关函数与功率谱密度
maxlag = 50;
[acf, lags] = xcorr(n_real, maxlag, 'unbiased');

[psd, f] = pwelch(n_real, hann(4096), 2048, 4096, fs, 'centered');

figure('Position', [100 100 980 400]);
subplot(1,2,1);
stem(lags, acf, 'b', 'filled', 'MarkerSize', 3);
grid on;
xlabel('时延 (样本)'); ylabel('R_n(\tau)');
title('自相关函数 (理论: \sigma^2\delta(\tau))');
subplot(1,2,2);
plot(f/1e3, 10*log10(psd), 'b');
hold on; yline(10*log10(sigma2/fs), 'r--', 'LineWidth', 1.5);
grid on; ylim([-60 -40]);
xlabel('频率 (kHz)'); ylabel('PSD (dB/Hz)');
legend('Welch估计', '理论白谱 \sigma^2/f_s', 'Location', 'southeast');
title('功率谱密度');
exportgraphics(gcf, '../figures/fig1_3_噪声自相关与功率谱.png', 'Resolution', 150);

fprintf('内容1完成, 图片已保存到 figures/ 目录\n');
