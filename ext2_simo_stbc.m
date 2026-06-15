%% 扩展题目2: 多天线系统性能仿真分析
% (1) SIMO QPSK: 单径瑞利衰落, 接收分集数 Nr=1/2/4,
%     合并方式: 最大比合并(MRC)/等增益合并(EGC)/选择合并(SC)
% (2) 4发1收STBC BPSK: 实正交满速率空时分组码, 与单发单收BPSK比较
clear; close all; rng(2026, 'twister');

%% 1. SIMO QPSK 不同分集数、不同合并方式
EbN0  = 0:2:30;
NrSet = [1 2 4];
combs = {'MRC', 'EGC', 'SC'};
minErr = 300; maxBits = 4e6; blk = 5e4;

berSIMO = zeros(numel(NrSet), numel(combs), numel(EbN0));
for in = 1:numel(NrSet)
    Nr = NrSet(in);
    for ie = 1:numel(EbN0)
        EsN0 = EbN0(ie) + 10*log10(2);
        N0   = 10^(-EsN0/10);
        nE = zeros(1,3); nB = 0;
        while min(nE) < minErr && nB < maxBits
            bits = randi([0 1], blk*2, 1);
            x = pskmod(bits, 4, pi/4, 'gray', 'InputType', 'bit');  % 1 x blk
            h = sqrt(0.5)*(randn(blk,Nr) + 1j*randn(blk,Nr));       % 单径瑞利
            n = sqrt(N0/2)*(randn(blk,Nr) + 1j*randn(blk,Nr));
            y = h .* x + n;                                          % 每支路接收
            % --- 三种合并 ---
            zMRC = sum(conj(h) .* y, 2);                             % 最大比
            zEGC = sum(exp(-1j*angle(h)) .* y, 2);                   % 等增益
            [~, sel] = max(abs(h), [], 2);                           % 选择合并
            idx  = sub2ind(size(h), (1:blk).', sel);
            zSC  = conj(h(idx)) .* y(idx);
            zAll = {zMRC, zEGC, zSC};
            for ic = 1:3
                bHat = pskdemod(zAll{ic}, 4, pi/4, 'gray', 'OutputType', 'bit');
                nE(ic) = nE(ic) + sum(bHat ~= bits);
            end
            nB = nB + numel(bits);
        end
        berSIMO(in, :, ie) = nE / nB;
        fprintf('Nr=%d Eb/N0=%2d dB: MRC=%.3e EGC=%.3e SC=%.3e\n', ...
                Nr, EbN0(ie), nE(1)/nB, nE(2)/nB, nE(3)/nB);
    end
end

% 图1: MRC不同分集数 + 理论
figure('Position', [100 100 720 480]);
sty = {'bo-', 'rs-', 'g^-'};
for in = 1:numel(NrSet)
    b = squeeze(berSIMO(in,1,:)); b(b==0) = nan;
    semilogy(EbN0, b, sty{in}, 'LineWidth', 1.4); hold on;
end
for in = 1:numel(NrSet)
    semilogy(EbN0, berfading(EbN0, 'psk', 4, NrSet(in)), 'k--', 'LineWidth', 0.9);
end
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('N_r=1 仿真', 'N_r=2 仿真', 'N_r=4 仿真', '理论(MRC)', 'Location', 'southwest');
title('SIMO QPSK最大比合并不同分集数的误比特率');
exportgraphics(gcf, '../figures/mimo_1_MRC分集数.png', 'Resolution', 150);

% 图2: 三种合并方式对比 (Nr=2 与 Nr=4)
figure('Position', [100 100 980 420]);
for sp = 1:2
    in = sp + 1;                       % Nr=2, Nr=4
    subplot(1,2,sp);
    for ic = 1:3
        b = squeeze(berSIMO(in,ic,:)); b(b==0) = nan;
        semilogy(EbN0, b, sty{ic}, 'LineWidth', 1.4); hold on;
    end
    semilogy(EbN0, squeeze(berfading(EbN0, 'psk', 4, NrSet(in))), 'k--');
    grid on; ylim([1e-6 1]);
    xlabel('E_b/N_0 (dB)'); ylabel('BER');
    legend('MRC', 'EGC', 'SC', 'MRC理论', 'Location', 'southwest');
    title(sprintf('N_r = %d', NrSet(in)));
end
sgtitle('SIMO QPSK三种合并方式误比特率对比');
exportgraphics(gcf, '../figures/mimo_2_合并方式对比.png', 'Resolution', 150);

%% 2. 4发1收实正交STBC BPSK vs 单发单收BPSK
% 实正交设计(4天线, 码率1): 行为时隙, 列为天线, 每符号块[x1 x2 x3 x4]
%   X = [ x1  x2  x3  x4
%        -x2  x1 -x4  x3
%        -x3  x4  x1 -x2
%        -x4 -x3  x2  x1 ] / 2     (总发射功率归一化)
% 接收: y = X h + n, 合并 s_hat = Re(A^H y), A^H A 实部 = ||h||^2 I
EbN0s = 0:2:30;
berSISO = zeros(size(EbN0s));
berSTBC = zeros(size(EbN0s));
nBlk = 5e4;                                     % 每批4符号块数
for ie = 1:numel(EbN0s)
    N0 = 10^(-EbN0s(ie)/10);                    % BPSK: Es = Eb (总发射能量)
    % ---- SISO BPSK ----
    nE = 0; nB = 0;
    while nE < minErr && nB < maxBits
        b  = randi([0 1], blk, 1);  x = 1 - 2*b;
        h  = sqrt(0.5)*(randn(blk,1)+1j*randn(blk,1));
        y  = h.*x + sqrt(N0/2)*(randn(blk,1)+1j*randn(blk,1));
        z  = real(conj(h).*y);
        nE = nE + sum((z<0) ~= b); nB = nB + blk;
    end
    berSISO(ie) = nE/nB;
    % ---- 4x1 STBC ----
    nE = 0; nB = 0;
    while nE < minErr && nB < maxBits
        b  = randi([0 1], nBlk, 4);  s = 1 - 2*b;          % 每行一个4符号块
        h  = sqrt(0.5)*(randn(nBlk,4)+1j*randn(nBlk,4));   % 块内准静态
        n  = sqrt(N0/2)*(randn(nBlk,4)+1j*randn(nBlk,4));
        % 4个时隙的接收信号(向量化): y_t = sum_i X(t,i) h_i / 2 + n_t
        y1 = ( s(:,1).*h(:,1) + s(:,2).*h(:,2) + s(:,3).*h(:,3) + s(:,4).*h(:,4))/2 + n(:,1);
        y2 = (-s(:,2).*h(:,1) + s(:,1).*h(:,2) - s(:,4).*h(:,3) + s(:,3).*h(:,4))/2 + n(:,2);
        y3 = (-s(:,3).*h(:,1) + s(:,4).*h(:,2) + s(:,1).*h(:,3) - s(:,2).*h(:,4))/2 + n(:,3);
        y4 = (-s(:,4).*h(:,1) - s(:,3).*h(:,2) + s(:,2).*h(:,3) + s(:,1).*h(:,4))/2 + n(:,4);
        % 正交合并 (A^H y 的实部)
        z1 = real( conj(h(:,1)).*y1 + conj(h(:,2)).*y2 + conj(h(:,3)).*y3 + conj(h(:,4)).*y4);
        z2 = real( conj(h(:,2)).*y1 - conj(h(:,1)).*y2 - conj(h(:,4)).*y3 + conj(h(:,3)).*y4);
        z3 = real( conj(h(:,3)).*y1 + conj(h(:,4)).*y2 - conj(h(:,1)).*y3 - conj(h(:,2)).*y4);
        z4 = real( conj(h(:,4)).*y1 - conj(h(:,3)).*y2 + conj(h(:,2)).*y3 - conj(h(:,1)).*y4);
        bHat = ([z1 z2 z3 z4] < 0);
        nE = nE + sum(bHat(:) ~= b(:)); nB = nB + numel(b);
    end
    berSTBC(ie) = nE/nB;
    fprintf('Eb/N0=%2d dB: SISO=%.3e, STBC4x1=%.3e\n', EbN0s(ie), berSISO(ie), berSTBC(ie));
end

figure('Position', [100 100 720 480]);
b1 = berSISO; b1(b1==0)=nan; b2 = berSTBC; b2(b2==0)=nan;
semilogy(EbN0s, b1, 'bo-', 'LineWidth', 1.4); hold on;
semilogy(EbN0s, berfading(EbN0s, 'psk', 2, 1), 'b--');
semilogy(EbN0s, b2, 'rs-', 'LineWidth', 1.4);
semilogy(EbN0s, berfading(EbN0s, 'psk', 2, 4), 'k-.', 'LineWidth', 1.0);
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('单发单收 仿真', '单发单收 理论', 'STBC 4\times1 仿真', ...
       '1\times4 MRC理论(参考)', 'Location', 'southwest');
title('4发1收实正交STBC BPSK与单发单收的误比特率比较');
exportgraphics(gcf, '../figures/mimo_3_STBC对比.png', 'Resolution', 150);

save('../figures/ber_mimo.mat', 'EbN0', 'NrSet', 'berSIMO', 'EbN0s', 'berSISO', 'berSTBC');
fprintf('扩展题目2(多天线)仿真完成\n');
