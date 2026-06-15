function [nErr, nBit, dbg] = ofdm_frame(EbN0dB, K, Ncp, doSync, P)
%OFDM_FRAME 仿真一帧OFDM/QPSK传输, 返回误比特数
%   EbN0dB : 信息比特信噪比(dB)
%   K      : 莱斯因子(线性), 定义为直射功率与全部散射功率之比
%            K=0 为纯瑞利多径, K=inf 为纯直射(AWGN式)
%   Ncp    : 循环前缀长度(样本)
%   doSync : true =随机定时偏移+CFO, 接收端两级同步(CP相关粗同步+
%            基于信道冲激响应的细同步)+LS信道估计
%            false=理想同步+理想信道状态信息
%   P      : 参数结构体(N, Nsym, tapDelay, tapPow, timBias)
%   dbg    : 调试信息(定时度量、信道响应等), 供演示绘图

N    = P.N;                          % 子载波数
Nsym = P.Nsym;                       % 每帧OFDM符号数(首符号为导频)
L    = N + Ncp;                      % 含CP符号长度
k    = 2;                            % QPSK每符号比特
EsN0 = EbN0dB + 10*log10(k);         % 子载波符号信噪比
N0   = 10^(-EsN0/10);

%% 发送端: 比特 -> QPSK -> IDFT -> 加CP
bitsTx = randi([0 1], N*k, Nsym);
X = pskmod(bitsTx(:), 4, pi/4, 'gray', 'InputType', 'bit');
X = reshape(X, N, Nsym);             % 频域符号, 第1列为导频(接收端已知)
x = ifft(X, N) * sqrt(N);            % IDFT, 保持时域平均功率为1
xt = [x(end-Ncp+1:end, :); x];       % 插入循环前缀
s  = xt(:);                          % 串行化

%% 莱斯多径信道(帧内准静态): 直射分量在首径, 总功率归一化为1
h  = zeros(max(P.tapDelay)+1, 1);
cn = sqrt(1/2) * (randn(numel(P.tapPow),1) + 1j*randn(numel(P.tapPow),1));
if isinf(K)
    h(1) = exp(1j*2*pi*rand);                       % 仅直射
else
    for ii = 1:numel(P.tapPow)
        h(P.tapDelay(ii)+1) = sqrt(P.tapPow(ii)/(K+1)) * cn(ii);
    end
    h(1) = h(1) + sqrt(K/(K+1)) * exp(1j*2*pi*rand);% 直射分量
end

if doSync
    off = randi([20 80]);            % 随机定时偏移
    eps_cfo = 0.2 * (2*rand - 1);    % 归一化载波频偏(子载波间隔为单位)
else
    off = 0; eps_cfo = 0;
end

tx = [zeros(off,1); s; zeros(N,1)];                    % 尾部留余量
rx = filter(h, 1, tx);                                 % 多径卷积
n  = (0:numel(rx)-1).';
rx = rx .* exp(1j*2*pi*eps_cfo*n/N);                   % 引入CFO
rx = rx + sqrt(N0/2)*(randn(size(rx)) + 1j*randn(size(rx)));

%% 接收端同步
dbg = struct();
if doSync
    % --- 第一级: CP相关粗定时 + CFO估计 ---
    dMax = 120;
    Pm = zeros(dMax+1, 1);
    for d = 0:dMax
        acc = 0;
        for m = 0:Nsym-1
            i1 = d + m*L + (1:Ncp).';
            if i1(end)+N > numel(rx), break; end
            acc = acc + sum(rx(i1) .* conj(rx(i1+N)));
        end
        Pm(d+1) = acc;
    end
    Ps = movmean(abs(Pm), 9);                          % 平滑抑制噪声抖动
    [~, dHat] = max(Ps);
    epsHat = -angle(Pm(dHat)) / (2*pi);                % CP相位估计CFO
    rx = rx .* exp(-1j*2*pi*epsHat*n/N);               % CFO校正
    start0 = max(dHat - 1 - P.timBias, 0);             % 粗定时(回退留余量)

    % --- 第二级: 基于LS信道冲激响应的细定时 ---
    Y1 = fft(rx(start0 + Ncp + (1:N))) / sqrt(N);
    hImp = ifft(Y1 ./ X(:,1));                         % 含循环移位的冲激响应
    e  = abs(hImp).^2;
    ew = movsum([e; e(1:Ncp)], [0 Ncp-1]);             % 长Ncp能量窗
    [~, t0] = max(ew(1:N));  t0 = t0 - 1;
    if t0 > N/2, t0 = t0 - N; end                      % 圆周折回负偏移
    start = max(start0 + t0 - 4, 0);                   % 首径前留4样本余量

    dbg.Pm = abs(Pm); dbg.off = off; dbg.dHat = dHat-1;
    dbg.epsTrue = eps_cfo; dbg.epsHat = epsHat;
    dbg.start0 = start0; dbg.start = start; dbg.t0 = t0;
else
    start = 0;
end

%% 解调: 去CP -> DFT -> 信道估计 -> 单抽头均衡
Y = zeros(N, Nsym);
for m = 0:Nsym-1
    seg = rx(start + m*L + Ncp + (1:N));
    Y(:, m+1) = fft(seg) / sqrt(N);
end

if doSync
    Hhat = Y(:,1) ./ X(:,1);                           % LS导频信道估计
else
    Hhat = fft(h, N);                                  % 理想CSI
end
dbg.Htrue = fft(h, N); dbg.Hhat = Hhat;

Xhat = Y(:, 2:end) ./ Hhat;                            % 迫零均衡(数据符号)
bitsRx = pskdemod(Xhat(:), 4, pi/4, 'gray', 'OutputType', 'bit');

ref  = bitsTx(:, 2:end);
nErr = sum(bitsRx ~= ref(:));
nBit = numel(ref);
end
