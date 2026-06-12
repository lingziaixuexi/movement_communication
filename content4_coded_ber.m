%% 内容4: 加入信道编码(RSC卷积码/LDPC码)与交织后的链路性能仿真
% 编码方案:
%   RSC : 码率1/2递归系统卷积码 (1, 15/13)_8, 约束长度4,
%         软判决Viterbi译码(vitdec 'unquant', 输入比特LLR)
%   LDPC: DVB-S.2标准LDPC码, 码率1/2, 码长64800, 和积译码25次迭代
% 交织: 随机交织器(randintrlv), 作用于编码后比特流
% 信道: 高斯信道 / 瑞利快衰落信道(理想CSI), 解调输出最大对数近似LLR
% 与内容3未编码系统在相同Eb/N0(按信息比特能量)下比较误码性能。
clear; close all; rng(2026, 'twister');

%% 公共参数
msgLen  = 32400;                       % 每帧信息比特数
R       = 1/2;                         % 码率
trellis = poly2trellis(4, [13 15], 13);% RSC (1, 15/13)_8
tblen   = 20;                          % Viterbi回溯深度
ldpcSeed= 9527;                        % 交织器种子
maxIter = 25;                          % LDPC最大迭代次数

Hldpc  = dvbs2ldpc(R);
cfgEnc = ldpcEncoderConfig(Hldpc);
cfgDec = ldpcDecoderConfig(Hldpc);

mods = struct('name', {'16QAM', '8PSK'}, 'M', {16, 8});
constels = { ...
    qammod((0:15).', 16, 'gray', 'UnitAveragePower', true), ...
    pskmod((0:7).', 8, pi/8, 'gray')};

EbN0_awgn = 0:1:9;                     % dB (编码系统)
EbN0_ray  = 0:2:20;                    % dB
minErr    = 100;                       % 每点最少错误比特
maxFrame  = 12;                        % 每点最大帧数

results = struct();

%% 主仿真循环: {16QAM,8PSK} x {AWGN,Rayleigh} x {RSC,LDPC}
for im = 1:2
    M = mods(im).M; k = log2(M);
    constel = constels{im};
    for ch = ["awgn", "ray"]
        if ch == "awgn", ebList = EbN0_awgn; else, ebList = EbN0_ray; end
        berRSC  = zeros(size(ebList));
        berLDPC = zeros(size(ebList));
        for ie = 1:numel(ebList)
            EsN0 = ebList(ie) + 10*log10(k*R);   % 信息比特能量折算
            N0   = 10^(-EsN0/10);
            errR = 0; errL = 0; nBit = 0;
            for fr = 1:maxFrame
                msg = randi([0 1], msgLen, 1);

                % ---- 编码 ----
                cwR = convenc(msg, trellis);          % RSC, 64800x1
                cwL = ldpcEncode(msg, cfgEnc);        % LDPC, 64800x1

                % ---- 交织 ----
                txR = randintrlv(cwR, ldpcSeed);
                txL = randintrlv(double(cwL), ldpcSeed);

                % ---- 调制 ----
                if M == 16
                    sR = qammod(txR, M, 'gray', 'InputType','bit', 'UnitAveragePower',true);
                    sL = qammod(txL, M, 'gray', 'InputType','bit', 'UnitAveragePower',true);
                else
                    sR = pskmod(txR, M, pi/M, 'gray', 'InputType','bit');
                    sL = pskmod(txL, M, pi/M, 'gray', 'InputType','bit');
                end
                Ns = numel(sR);

                % ---- 信道 ----
                if ch == "awgn"
                    h = ones(Ns, 1);
                else
                    h = sqrt(0.5) * (randn(Ns,1) + 1j*randn(Ns,1));  % 快衰落
                end
                n  = sqrt(N0/2) * (randn(Ns,1) + 1j*randn(Ns,1));
                yR = h .* sR + n;
                yL = h .* sL + n;

                % ---- 软解调 + 解交织 ----
                llrR = randdeintrlv(soft_llr(yR, h, constel, N0), ldpcSeed);
                llrL = randdeintrlv(soft_llr(yL, h, constel, N0), ldpcSeed);

                % ---- 译码 ----
                msgR = vitdec(llrR, trellis, tblen, 'trunc', 'unquant');
                msgL = double(ldpcDecode(llrL, cfgDec, maxIter));

                errR = errR + sum(msgR ~= msg);
                errL = errL + sum(msgL ~= msg);
                nBit = nBit + msgLen;
                if errR >= minErr && errL >= minErr && fr >= 3
                    break;
                end
            end
            berRSC(ie)  = errR / nBit;
            berLDPC(ie) = errL / nBit;
            fprintf('%s-%s Eb/N0=%2d dB: RSC BER=%.3e, LDPC BER=%.3e (%d bits)\n', ...
                    mods(im).name, ch, ebList(ie), berRSC(ie), berLDPC(ie), nBit);
        end
        results.(sprintf('m%d_%s_rsc',  M, ch)) = berRSC;
        results.(sprintf('m%d_%s_ldpc', M, ch)) = berLDPC;
    end
end

%% 绘图: 与未编码理论值对比
figCfg = {                       % {调制idx, 信道, Eb/N0列表, 图号}
    1, "awgn", EbN0_awgn, 1;
    1, "ray",  EbN0_ray,  2;
    2, "awgn", EbN0_awgn, 3;
    2, "ray",  EbN0_ray,  4};
chName = struct('awgn', '高斯信道', 'ray', '瑞利衰落信道');

for ic = 1:size(figCfg, 1)
    im = figCfg{ic,1}; ch = figCfg{ic,2}; eb = figCfg{ic,3}; figNo = figCfg{ic,4};
    M  = mods(im).M;
    if ch == "awgn"
        if M == 16, thUn = berawgn(eb, 'qam', 16);
        else,       thUn = berawgn(eb, 'psk', 8, 'nondiff'); end
    else
        if M == 16, thUn = berfading(eb, 'qam', 16, 1);
        else,       thUn = berfading(eb, 'psk', 8, 1); end
    end
    bR = results.(sprintf('m%d_%s_rsc',  M, ch));
    bL = results.(sprintf('m%d_%s_ldpc', M, ch));
    bR(bR == 0) = nan; bL(bL == 0) = nan;

    figure('Position', [100 100 720 480]);
    semilogy(eb, thUn, 'k--', 'LineWidth', 1.2); hold on;
    semilogy(eb, bR, 'bo-', 'LineWidth', 1.4);
    semilogy(eb, bL, 'rs-', 'LineWidth', 1.4);
    grid on; ylim([1e-6 1]);
    xlabel('E_b/N_0 (dB)'); ylabel('BER');
    legend('未编码(理论)', 'RSC(1,15/13)+交织', 'LDPC(1/2)+交织', ...
           'Location', 'southwest');
    title(sprintf('%s %s 编码前后误码性能', mods(im).name, chName.(ch)));
    exportgraphics(gcf, sprintf('../figures/fig4_%d_%s_%s编码对比.png', ...
                   figNo, mods(im).name, chName.(ch)), 'Resolution', 150);
end

%% 交织器作用演示: 时间相关瑞利衰落(突发错误)下 RSC+16QAM 有无交织对比
fdTs    = 5e-4;                        % 归一化多普勒, 衰落相干时间约846符号
ebList  = 0:2:20;
M = 16; k = 4; constel = constels{1};
berNoInt = zeros(size(ebList));
berInt   = zeros(size(ebList));
for ie = 1:numel(ebList)
    EsN0 = ebList(ie) + 10*log10(k*R);
    N0   = 10^(-EsN0/10);
    eN = 0; eI = 0; nBit = 0;
    for fr = 1:maxFrame
        msg = randi([0 1], msgLen, 1);
        cw  = convenc(msg, trellis);
        txN = cw;                             % 无交织
        txI = randintrlv(cw, ldpcSeed);       % 有交织
        sN  = qammod(txN, M, 'gray', 'InputType','bit', 'UnitAveragePower',true);
        sI  = qammod(txI, M, 'gray', 'InputType','bit', 'UnitAveragePower',true);
        Ns  = numel(sN);
        h   = jakes_rayleigh(Ns, 1, fdTs, 32);     % 时间相关衰落
        n   = sqrt(N0/2) * (randn(Ns,1) + 1j*randn(Ns,1));
        llrN = soft_llr(h.*sN + n, h, constel, N0);
        llrI = randdeintrlv(soft_llr(h.*sI + n, h, constel, N0), ldpcSeed);
        mN  = vitdec(llrN, trellis, tblen, 'trunc', 'unquant');
        mI  = vitdec(llrI, trellis, tblen, 'trunc', 'unquant');
        eN  = eN + sum(mN ~= msg);
        eI  = eI + sum(mI ~= msg);
        nBit = nBit + msgLen;
        if eN >= minErr && eI >= minErr && fr >= 3, break; end
    end
    berNoInt(ie) = eN / nBit;
    berInt(ie)   = eI / nBit;
    fprintf('交织对比 Eb/N0=%2d dB: 无交织 BER=%.3e, 有交织 BER=%.3e\n', ...
            ebList(ie), berNoInt(ie), berInt(ie));
end

berNoInt(berNoInt == 0) = nan; berInt(berInt == 0) = nan;
figure('Position', [100 100 720 480]);
semilogy(ebList, berfading(ebList, 'qam', 16, 1), 'k--', 'LineWidth', 1.2); hold on;
semilogy(ebList, berNoInt, 'mv-', 'LineWidth', 1.4);
semilogy(ebList, berInt, 'bo-', 'LineWidth', 1.4);
grid on; ylim([1e-6 1]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('未编码(理论)', 'RSC无交织', 'RSC+随机交织', 'Location', 'southwest');
title(sprintf('时间相关瑞利信道下交织的作用 (16QAM, f_dT_s=%.0e)', fdTs));
exportgraphics(gcf, '../figures/fig4_5_交织作用对比.png', 'Resolution', 150);

save('../figures/ber_coded.mat', 'results', 'EbN0_awgn', 'EbN0_ray', ...
     'ebList', 'berNoInt', 'berInt');
fprintf('内容4完成, 图片已保存到 figures/ 目录\n');
