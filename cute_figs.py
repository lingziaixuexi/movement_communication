# -*- coding: utf-8 -*-
"""可爱手绘风(xkcd/excalidraw)原理插图: OFDM子载波正交性、莱斯分布、分集原理、循环前缀"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = "../figures/"
CUTE = ["Hannotate SC", "HanziPen SC", "PingFang SC"]

def cute_style():
    ctx = plt.xkcd(scale=1.1, length=120, randomness=1.6)
    plt.rcParams["font.family"] = CUTE
    plt.rcParams["axes.unicode_minus"] = False
    return ctx

# ---------- 1. OFDM 子载波正交性 ----------
with cute_style():
    fig, ax = plt.subplots(figsize=(8.6, 4.4), dpi=150)
    f = np.linspace(-1.5, 6.5, 1200)
    colors = ["#e07a5f", "#3d7ea6", "#81b29a", "#f2a541", "#9b5de5"]
    for k in range(5):
        y = np.sinc(f - k)
        ax.plot(f, y, lw=2, color=colors[k])
    ax.axhline(0, color="gray", lw=1)
    for k in range(5):
        ax.plot(k, 1, "o", ms=7, color=colors[k])
        ax.plot(k, 0, "x", ms=8, mew=2, color="dimgray")
    ax.annotate("某个子载波取峰值时,\n其余子载波恰好过零!",
                xy=(2, 1.0), xytext=(3.9, 1.05),
                arrowprops=dict(arrowstyle="->", color="black"),
                fontsize=12)
    ax.annotate("频谱可以重叠\n却互不干扰",
                xy=(0.5, 0.62), xytext=(-1.4, 0.88),
                arrowprops=dict(arrowstyle="->", color="black"),
                fontsize=12)
    ax.set_xlabel("频率 (以子载波间隔 $\\Delta f$ 为单位)", fontsize=13)
    ax.set_ylabel("幅度", fontsize=13)
    ax.set_title("OFDM 子载波的正交性", fontsize=15)
    ax.set_ylim(-0.35, 1.5)
    ax.set_yticks([0, 0.5, 1])
    fig.tight_layout()
    fig.savefig(OUT + "cute_ofdm_subcarriers.png")
    plt.close(fig)

# ---------- 2. 莱斯分布 PDF 随 K 变化 ----------
with cute_style():
    fig, ax = plt.subplots(figsize=(8.6, 4.4), dpi=150)
    r = np.linspace(0, 3.2, 600)
    from numpy import i0
    for K, c in zip([0, 1, 3, 10], ["#3d7ea6", "#81b29a", "#f2a541", "#e07a5f"]):
        # 总功率归一化: 2*sigma^2*(1+K)=1
        s2 = 1 / (2 * (1 + K))
        A = np.sqrt(K / (K + 1))
        pdf = (r / s2) * np.exp(-(r**2 + A**2) / (2 * s2)) * i0(r * A / s2)
        ax.plot(r, pdf, lw=2, color=c, label=f"K={K}")
    ax.set_ylim(0, 2.45)
    ax.annotate("K=0 就是瑞利分布\n(没有直射径, 衰落最狠)",
                xy=(0.5, 0.83), xytext=(0.05, 1.6),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=12)
    ax.annotate("K 越大越集中,\n越接近不衰落的AWGN",
                xy=(1.07, 1.9), xytext=(1.75, 1.45),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=12)
    ax.set_xlabel("包络 r", fontsize=13)
    ax.set_ylabel("概率密度 f(r)", fontsize=13)
    ax.set_title("莱斯分布: 直射径越强(K越大), 信道越稳", fontsize=15)
    ax.legend(fontsize=11, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT + "cute_rician_pdf.png")
    plt.close(fig)

# ---------- 3. 分集原理 ----------
with cute_style():
    rng = np.random.default_rng(7)
    fig, ax = plt.subplots(figsize=(8.6, 4.4), dpi=150)
    t = np.linspace(0, 1, 800)
    env = []
    for i, c in zip(range(2), ["#3d7ea6", "#81b29a"]):
        g = (rng.standard_normal(9) + 1j * rng.standard_normal(9)) / np.sqrt(2)
        fd = np.arange(1, 10) * 2.2
        ph = rng.uniform(0, 2 * np.pi, 9)
        h = np.sum([gk * np.exp(1j * (2 * np.pi * fk * t + pk))
                    for gk, fk, pk in zip(g, fd, ph)], axis=0) / 3
        e = 20 * np.log10(np.abs(h) + 1e-3)
        env.append(e)
        ax.plot(t, e, lw=1.6, color=c, alpha=0.85, label=f"支路 {i+1}")
    best = np.maximum(env[0], env[1])
    ax.plot(t, best, lw=2.8, color="#e07a5f", label="选大的那条!")
    deep = np.argmin(env[0])
    ax.annotate("支路1掉进深衰落坑里",
                xy=(t[deep], max(env[0][deep], -38)), xytext=(0.03, -27),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=12)
    ax.annotate("但支路2没掉坑,\n合并后稳稳的~",
                xy=(t[deep], best[deep]), xytext=(0.6, 8),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=12)
    ax.set_xlabel("时间", fontsize=13)
    ax.set_ylabel("包络 (dB)", fontsize=13)
    ax.set_title("分集的思想: 多条独立衰落支路同时掉坑的概率很小", fontsize=14)
    ax.set_ylim(-40, 18)
    ax.legend(fontsize=11, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT + "cute_diversity.png")
    plt.close(fig)

# ---------- 4. 循环前缀对抗多径 ----------
with cute_style():
    fig, ax = plt.subplots(figsize=(8.6, 4.0), dpi=150)
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.4); ax.axis("off")

    def sym(x, y, wcp=0.9, wu=3.4, color="#3d7ea6", label=True):
        ax.add_patch(Rectangle((x, y), wcp, 0.8, facecolor="#f2a541",
                               edgecolor="k", lw=1.5, hatch="//"))
        ax.add_patch(Rectangle((x + wcp, y), wu, 0.8, facecolor=color,
                               edgecolor="k", lw=1.5, alpha=0.55))
        if label:
            ax.text(x + wcp / 2, y + 0.4, "CP", ha="center", va="center", fontsize=11)
            ax.text(x + wcp + wu / 2, y + 0.4, "OFDM符号正文", ha="center",
                    va="center", fontsize=12)

    # 直射径
    sym(1.2, 3.9); sym(5.5, 3.9)
    ax.text(0.15, 4.3, "直射径", fontsize=12, va="center")
    # 时延径
    sym(1.9, 2.6); sym(6.2, 2.6)
    ax.text(0.15, 3.0, "时延径", fontsize=12, va="center")
    # 复制示意: 正文尾部 -> CP
    arr = FancyArrowPatch((5.4, 4.9), (1.6, 4.9), arrowstyle="->",
                          mutation_scale=14, color="#e07a5f",
                          connectionstyle="arc3,rad=0.25", lw=1.8)
    ax.add_patch(arr)
    ax.text(3.5, 5.15, "把正文尾巴复制到开头 = 循环前缀", fontsize=12,
            ha="center", color="#e07a5f")
    # FFT 窗
    ax.add_patch(Rectangle((2.1, 2.45), 3.4, 2.5, fill=False,
                           edgecolor="#9b5de5", lw=2.2, linestyle="--"))
    ax.text(3.8, 2.1, "FFT窗里只有完整的循环内容,\n上一符号的回声只打在CP上→ISI被吸收!",
            fontsize=12, ha="center", va="top")
    ax.set_title("循环前缀: 给多径回声准备的\"缓冲垫\"", fontsize=15)
    fig.tight_layout()
    fig.savefig(OUT + "cute_cp.png")
    plt.close(fig)

print("4张可爱风插图已生成")
