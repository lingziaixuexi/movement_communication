%% 内容3: 16QAM与8PSK调制解调及高斯/瑞利信道误码率仿真
% (1) 观察16QAM/8PSK信号经过高斯信道、瑞利衰落信道前后的星座图;
% (2) 统计两种调制在两种信道下的误比特率, 并与理论曲线对比。
% 瑞利信道下接收端假设理想信道估计, 采用相干检测(迫零均衡 y/h)。
clear; close all; rng(2026, 'twister');

mods = struct( ...
    'name', {'16QAM', '8PSK'}, ...
    'M',    {16, 8});

%% 1. 星座图观察 (Es/N0 = 18 dB)
EsN0_demo = 18;                       % dB
Nsym_demo = 3000;
fdTs      = 5e-4;                     % 星座演示用时间相关衰落

for im = 1:2
    M = mods(im).M; k = log2(M);
    bits = randi([0 1], Nsym_demo*k, 1);
    if M == 16
        s = qammod(bits, M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
    else
        s = pskmod(bits, M, pi/M, 'gray', 'InputType', 'bit');
    end
    N0 = 10^(-EsN0_demo/10);          % Es = 1
    n  = sqrt(N0/2) * (randn(Nsym_demo,1) + 1j*randn(Nsym_demo,1));
    h  = jakes_rayleigh(Nsym_demo, 1, fdTs, 32);   % 归一化采样率

    y_awgn = s + n;                   % 高斯信道
    y_ray  = h .* s + n;              % 瑞利衰落信道
    y_eq   = y_ray ./ h;              % 理想信道估计后迫零均衡

    figure('Position', [80 80 900 760]);
    subplot(2,2,1);
    plot(real(s), imag(s), 'b.'); grid on; axis equal; axis([-2 2 -2 2]);
    title('发送星座'); xlabel('I'); ylabel('Q');
    subplot(2,2,2);
    plot(real(y_awgn), imag(y_awgn), 'b.'); grid on; axis equal; axis([-2 2 -2 2]);
    title(sprintf('经高斯信道 (E_s/N_0=%d dB)', EsN0_demo)); xlabel('I'); ylabel('Q');
    subplot(2,2,3);
    plot(real(y_ray), imag(y_ray), 'b.'); grid on; axis equal; axis([-2 2 -2 2]);
    title('经瑞利衰落信道(均衡前)'); xlabel('I'); ylabel('Q');
    subplot(2,2,4);
    plot(real(y_eq), imag(y_eq), 'b.'); grid on; axis equal; axis([-2 2 -2 2]);
    title('瑞利信道理想均衡后'); xlabel('I'); ylabel('Q');
    sgtitle(sprintf('%s 星座图', mods(im).name));
    exportgraphics(gcf, sprintf('../figures/fig3_%d_%s星座图.png', im, mods(im).name), ...
                   'Resolution', 150);
end

%% 2. 误比特率仿真
EbN0_awgn = 0:2:16;                   % dB
EbN0_ray  = 0:4:36;                   % dB
minErr  = 300;                        % 每点最少错误比特数
maxBits = 4e6;                        % 每点最大仿真比特数
blkSym  = 5e4;                        % 每批符号数

ber = struct();
for im = 1:2
    M = mods(im).M; k = log2(M);
    for ch = ["awgn", "ray"]
        if ch == "awgn", ebList = EbN0_awgn; else, ebList = EbN0_ray; end
        berSim = zeros(size(ebList));
        for ie = 1:numel(ebList)
            EsN0 = ebList(ie) + 10*log10(k);
            N0   = 10^(-EsN0/10);
            nErr = 0; nBit = 0;
            while nErr < minErr && nBit < maxBits
                bits = randi([0 1], blkSym*k, 1);
                if M == 16
                    s = qammod(bits, M, 'gray', 'InputType','bit', 'UnitAveragePower',true);
                else
                    s = pskmod(bits, M, pi/M, 'gray', 'InputType','bit');
                end
                n = sqrt(N0/2) * (randn(blkSym,1) + 1j*randn(blkSym,1));
                if ch == "awgn"
                    y = s + n; heq = ones(blkSym,1);
                else
                    heq = sqrt(0.5) * (randn(blkSym,1) + 1j*randn(blkSym,1)); % 快衰落
                    y = heq .* s + n;
                end
                yEq = y ./ heq;
                if M == 16
                    bHat = qamdemod(yEq, M, 'gray', 'OutputType','bit', ...
                                    'UnitAveragePower', true);
                else
                    bHat = pskdemod(yEq, M, pi/M, 'gray', 'OutputType','bit');
                end
                nErr = nErr + sum(bHat ~= bits);
                nBit = nBit + numel(bits);
            end
            berSim(ie) = nErr / nBit;
            fprintf('%s-%s Eb/N0=%2d dB: BER = %.3e (%d bits)\n', ...
                    mods(im).name, ch, ebList(ie), berSim(ie), nBit);
        end
        ber.(sprintf('m%d_%s', M, ch)) = berSim;
    end
end

% 理论值
th16_awgn = berawgn(EbN0_awgn, 'qam', 16);
th8_awgn  = berawgn(EbN0_awgn, 'psk', 8, 'nondiff');
th16_ray  = berfading(EbN0_ray, 'qam', 16, 1);
th8_ray   = berfading(EbN0_ray, 'psk', 8, 1);

%% 3. BER曲线
figure('Position', [100 100 980 420]);
subplot(1,2,1);
semilogy(EbN0_awgn, ber.m16_awgn, 'bo-', 'LineWidth', 1.2); hold on;
semilogy(EbN0_awgn, th16_awgn, 'b--');
semilogy(EbN0_awgn, ber.m8_awgn, 'rs-', 'LineWidth', 1.2);
semilogy(EbN0_awgn, th8_awgn, 'r--');
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('16QAM 仿真', '16QAM 理论', '8PSK 仿真', '8PSK 理论', 'Location', 'southwest');
title('高斯信道');
subplot(1,2,2);
semilogy(EbN0_ray, ber.m16_ray, 'bo-', 'LineWidth', 1.2); hold on;
semilogy(EbN0_ray, th16_ray, 'b--');
semilogy(EbN0_ray, ber.m8_ray, 'rs-', 'LineWidth', 1.2);
semilogy(EbN0_ray, th8_ray, 'r--');
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('16QAM 仿真', '16QAM 理论', '8PSK 仿真', '8PSK 理论', 'Location', 'southwest');
title('瑞利衰落信道(理想CSI相干检测)');
exportgraphics(gcf, '../figures/fig3_3_未编码BER曲线.png', 'Resolution', 150);

save('../figures/ber_uncoded.mat', 'ber', 'EbN0_awgn', 'EbN0_ray');
fprintf('内容3完成, 图片已保存到 figures/ 目录\n');
