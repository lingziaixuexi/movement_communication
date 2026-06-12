function llr = soft_llr(y, h, constel, noiseVar)
%SOFT_LLR 任意星座的最大对数近似比特软信息(LLR)计算
%   llr = soft_llr(y, h, constel, noiseVar)
%   y        : N x 1 接收符号
%   h        : N x 1 信道增益(AWGN信道置 1), 假设接收端理想信道估计
%   constel  : M x 1 星座点, 第 m 个点对应整数 m-1 (MSB在前) 的比特组
%   noiseVar : 复噪声总方差 N0 (标量或 N x 1 向量)
%   输出 llr : (N*log2(M)) x 1, 约定 LLR = log[P(b=0)/P(b=1)], 正值判 0
%
%   LLR_b ≈ ( min_{s∈S1} |y - h*s|^2 - min_{s∈S0} |y - h*s|^2 ) / N0

M = numel(constel);
k = log2(M);
% M 个星座点对应的比特图样 (M x k, MSB 在前)
bmat = reshape(int2bit((0:M-1).', k), k, M).';

d2 = abs(y - h .* constel.').^2 ./ noiseVar;   % N x M 归一化欧氏距离平方
llrMat = zeros(numel(y), k);
for b = 1:k
    s0 = (bmat(:, b) == 0);
    llrMat(:, b) = min(d2(:, ~s0), [], 2) - min(d2(:, s0), [], 2);
end
llr = reshape(llrMat.', [], 1);   % 按符号串行化, 每符号 MSB 在前
end
