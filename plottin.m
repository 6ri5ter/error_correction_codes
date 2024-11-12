clc
clear
close all
% M = 2.^(1:9);
% SNR = 10.^((0:5:35)/10);
% m = 2:14;
% modul_name = {'PAM','PSK','QAM'};
%[data,data_rate] = parameter_sweep(m,M,modul_name,SNR);
data = load('one_sim_trta.mat');
for i = 1:3
    for j = [3,5,7]
        ecc_plot_1(data,{'block','d','M'},i,[1:9],j)
    end
end