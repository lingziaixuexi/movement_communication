%% 一键运行基本题目全部仿真 (内容1~4)
% 运行前确保当前目录为 code/, 图片输出到 ../figures/
tStart = tic;
content1_gaussian_noise;
content2_rayleigh_channel;
content3_qam_psk;
content4_coded_ber;
fprintf('\n全部仿真完成, 总耗时 %.1f 分钟\n', toc(tStart)/60);
