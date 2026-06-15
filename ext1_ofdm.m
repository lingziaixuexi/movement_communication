%% 扩展题目1: OFDM系统性能仿真
% 设计依据: 最大时延扩展 tau_max = 4 us, 莱斯3径信道, 莱斯因子K=3
%   采样率 fs = 5 MHz (Ts = 0.2 us), 子载波数 N = 512
%   子载波间隔 df = fs/N = 9.766 kHz << 相干带宽 Bc ~ 143 kHz
%   有用符号周期 Tu = 1/df = 102.4 us >> tau_max
%   CP长度 Ncp = 32 样本 (6.4 us > tau_max = 4 us), 开销 6.25%
%   总符号周期 Tsym = (512+32)*0.2us = 108.8 us
% 链路: QPSK -> IDFT -> 加CP -> 莱斯3径信道 -> CP相关定时同步/CFO估计
%       -> DFT -> LS信道估计 -> 迫零均衡 -> 判决
clear; close all; rng(2026, 'twister');

P.N        = 512;
P.Nsym     = 10;                          % 1导频 + 9数据
P.fs       = 5e6;
P.tapDelay = [0 8 20];                    % 0, 1.6, 4 us @5MHz
tapPowdB   = [0 -3 -6];
P.tapPow   = 10.^(tapPowdB/10) / sum(10.^(tapPowdB/10));
P.timBias  = 11;                          % 定时回退量(样本), 利用CP余量
Ncp        = 32;

%% 1. 同步与信道估计演示 (K=3, Eb/N0=15 dB)
rng(7);
[~, ~, dbg] = ofdm_frame(15, 3, Ncp, true, P);
fprintf('定时: 真实偏移=%d, 估计=%d (提前%d样本使用)\n', dbg.off, dbg.dHat, P.timBias);
fprintf('CFO : 真实=%.4f, 估计=%.4f (子载波间隔归一化)\n', dbg.epsTrue, dbg.epsHat);

figure('Position', [100 100 980 400]);
subplot(1,2,1);
plot(0:numel(dbg.Pm)-1, dbg.Pm/max(dbg.Pm), 'b', 'LineWidth', 1.2); hold on;
xline(dbg.off, 'r--', 'LineWidth', 1.2);
plot(dbg.dHat, 1, 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
xlabel('候选定时位置 d (样本)'); ylabel('归一化定时度量 |P(d)|');
legend('CP相关度量', '真实定时位置', '估计峰值', 'Location', 'southeast');
title('CP相关定时同步度量');
subplot(1,2,2);
fsc = (0:P.N-1).';
plot(fsc, 20*log10(abs(dbg.Htrue)), 'b', 'LineWidth', 1.2); hold on;
plot(fsc, 20*log10(abs(dbg.Hhat)), 'r:', 'LineWidth', 1.2);
grid on; ylim([-30 10]);
xlabel('子载波序号 k'); ylabel('|H(k)| (dB)');
legend('真实信道响应', 'LS估计', 'Location', 'southwest');
title('信道频率响应与LS估计 (E_b/N_0=15 dB)');
exportgraphics(gcf, '../figures/ofdm_1_同步与信道估计.png', 'Resolution', 150);

%% 2. 不同莱斯因子K下的BER (理想同步+理想CSI, 隔离K的影响)
rng(2026);
EbN0  = 0:2:30;
Klist = [0 3 10];
minErr = 300; maxFrm = 800; minFrm = 200;  % minFrm保证信道实现遍历性
berK = zeros(numel(Klist), numel(EbN0));
for ik = 1:numel(Klist)
    for ie = 1:numel(EbN0)
        nE = 0; nB = 0; fr = 0;
        while (nE < minErr || fr < minFrm) && fr < maxFrm
            [e, b] = ofdm_frame(EbN0(ie), Klist(ik), Ncp, false, P);
            nE = nE + e; nB = nB + b; fr = fr + 1;
        end
        berK(ik, ie) = nE / nB;
        fprintf('K=%2d Eb/N0=%2d dB: BER=%.3e (%d帧)\n', Klist(ik), EbN0(ie), berK(ik,ie), fr);
    end
end

figure('Position', [100 100 720 480]);
styles = {'bo-', 'rs-', 'g^-'};
for ik = 1:numel(Klist)
    b = berK(ik,:); b(b==0) = nan;
    if Klist(ik) == 10, b(EbN0 > 20) = nan; end   % 高信噪比点错误事件不足
    semilogy(EbN0, b, styles{ik}, 'LineWidth', 1.4); hold on;
end
semilogy(EbN0, berawgn(EbN0, 'psk', 4, 'nondiff'), 'k--', 'LineWidth', 1.2);
semilogy(EbN0, berfading(EbN0, 'psk', 4, 1), 'm-.', 'LineWidth', 1.0);
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('K=0 (瑞利)', 'K=3', 'K=10', 'AWGN理论', '平坦瑞利理论', 'Location', 'southwest');
title('OFDM/QPSK不同莱斯因子下的误比特率');
exportgraphics(gcf, '../figures/ofdm_2_不同K因子BER.png', 'Resolution', 150);

%% 3. 实际同步+LS估计 vs 理想同步+理想CSI (K=3)
rng(2027);
minErr2 = 600; minFrm2 = 300; maxFrm2 = 1500;
berIdeal = zeros(size(EbN0));  berReal = zeros(size(EbN0));
for ie = 1:numel(EbN0)
    cfgs = {false, true};
    out  = zeros(1,2);
    for ic = 1:2
        nE = 0; nB = 0; fr = 0;
        while (nE < minErr2 || fr < minFrm2) && fr < maxFrm2
            [e, b] = ofdm_frame(EbN0(ie), 3, Ncp, cfgs{ic}, P);
            nE = nE + e; nB = nB + b; fr = fr + 1;
        end
        out(ic) = nE / nB;
    end
    berIdeal(ie) = out(1); berReal(ie) = out(2);
    fprintf('K=3 Eb/N0=%2d dB: 理想=%.3e, 同步+LS=%.3e\n', EbN0(ie), out(1), out(2));
end

figure('Position', [100 100 720 480]);
b1 = berIdeal; b1(b1==0)=nan; b2 = berReal; b2(b2==0)=nan;
semilogy(EbN0, b1, 'bo-', 'LineWidth', 1.4); hold on;
semilogy(EbN0, b2, 'rs-', 'LineWidth', 1.4);
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('理想同步+理想CSI', '两级同步+LS信道估计', 'Location', 'southwest');
title('OFDM/QPSK实际同步与信道估计的性能代价 (K=3)');
exportgraphics(gcf, '../figures/ofdm_3_同步与估计代价.png', 'Resolution', 150);

%% 3b. CP长度的影响 (K=0纯瑞利多径, 理想同步+理想CSI)
rng(2028);
cpList = [32 4 0];
berCP = zeros(numel(cpList), numel(EbN0));
for ic = 1:numel(cpList)
    for ie = 1:numel(EbN0)
        nE = 0; nB = 0; fr = 0;
        while (nE < minErr || fr < minFrm) && fr < maxFrm
            [e, b] = ofdm_frame(EbN0(ie), 0, cpList(ic), false, P);
            nE = nE + e; nB = nB + b; fr = fr + 1;
        end
        berCP(ic, ie) = nE / nB;
    end
    fprintf('Ncp=%2d 完成\n', cpList(ic));
end

figure('Position', [100 100 720 480]);
styles2 = {'bo-', 'mv-', 'rs-'};
lbl = {'N_{cp}=32 (6.4\mus > \tau_{max})', 'N_{cp}=4 (0.8\mus < \tau_{max})', '无CP'};
for ic = 1:numel(cpList)
    b = berCP(ic,:); b(b==0) = nan;
    semilogy(EbN0, b, styles2{ic}, 'LineWidth', 1.4); hold on;
end
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend(lbl, 'Location', 'southwest');
title('CP长度对OFDM/QPSK性能的影响 (K=0)');
exportgraphics(gcf, '../figures/ofdm_5_CP长度影响.png', 'Resolution', 150);

%% 4. 均衡前后星座图 (K=3, Eb/N0=15 dB)
rng(11);
NsymC = 6;  Pc = P; Pc.Nsym = NsymC;
bits = randi([0 1], P.N*2, NsymC);
X = reshape(pskmod(bits(:), 4, pi/4, 'gray', 'InputType','bit'), P.N, NsymC);
x = ifft(X, P.N)*sqrt(P.N); xt = [x(end-Ncp+1:end,:); x]; s = xt(:);
h = zeros(21,1);                          % K=3莱斯信道(同ofdm_frame约定)
cn = sqrt(1/2)*(randn(3,1)+1j*randn(3,1));
Kc = 3;
h([1 9 21]) = sqrt(P.tapPow/(Kc+1)).' .* cn;
h(1) = h(1) + sqrt(Kc/(Kc+1))*exp(1j*2*pi*rand);
N0 = 10^(-(15+10*log10(2))/10);
rx = filter(h,1,s) + sqrt(N0/2)*(randn(size(s))+1j*randn(size(s)));
Y = zeros(P.N, NsymC);
for m = 0:NsymC-1
    Y(:,m+1) = fft(rx(m*(P.N+Ncp)+Ncp+(1:P.N)))/sqrt(P.N);
end
H = fft(h, P.N);  Yeq = Y ./ H;
figure('Position', [100 100 900 420]);
subplot(1,2,1);
plot(real(Y(:)), imag(Y(:)), 'b.', 'MarkerSize', 4); grid on; axis equal;
xlabel('I'); ylabel('Q'); title('均衡前 (频率选择性衰落)');
subplot(1,2,2);
plot(real(Yeq(:)), imag(Yeq(:)), 'b.', 'MarkerSize', 4); grid on; axis equal; axis([-3 3 -3 3]);
xlabel('I'); ylabel('Q'); title('单抽头迫零均衡后');
sgtitle('OFDM/QPSK子载波星座 (K=3, E_b/N_0=15 dB)');
exportgraphics(gcf, '../figures/ofdm_4_均衡前后星座.png', 'Resolution', 150);

save('../figures/ber_ofdm.mat', 'EbN0', 'Klist', 'berK', 'berIdeal', 'berReal', 'berCP', 'cpList');
fprintf('扩展题目1(OFDM)仿真完成\n');
