%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% Demo for how to use the suite
clear
clc
close
rng('default')
[ri_opt,li_opt,rate_opt,rx,lx,rate] = li_ri_opt(16,6,100);
li = li_opt{1};
ri = ri_opt{1};
opt_rate = rate_opt;
ns = [100,1000,5000,10000,20000];
epsilons = 0.05:0.05:0.3;
block_err = zeros(length(ns),length(epsilons));
ber = zeros(length(ns),length(epsilons));
numofbits = 1e5;
for i = 1:length(ns)
    fprintf('\nn = %d\n',ns(i))
    [Rf,Lf,newrate] = Li_Ri_approx_opt_v1_2(ri,li,ns(i));
    rmax = find(Rf>0,1,'last');
    lmax = find(Lf>0,1,'last');
    Rf = Rf(1:rmax);
    Lf = Lf(1:lmax);
    H = poly2mat(Rf,Lf);
    H(:,~sum(H)) = [];
    n = size(H,2);
    blocks = ceil(numofbits/n);
    data = zeros(blocks,n);%encoded data
    for j = 1:length(epsilons)
        tic
        fprintf('e = %d\n',epsilons(j))
        channel_data = data;%simulating bec
        add_er = rand(blocks,n);
        channel_data(add_er<epsilon) = nan;
        corrected = bp_dec(channel_data,H,epsilons(j));
        %Calculating Performance
        iscor = zeros(size(data,1),1);
        corr_bits = iscor;
        for k = 1:size(data,1)
            corr_bits(k) = sum(corrected(k,:)==data(k,:));
            iscor(k) = corr_bits(k)==n;
        end
        block_err(i,j) = 1-sum(iscor)/length(iscor);
        ber(i,j) = 1-sum(corr_bits)/numel(data);
        toc
    end
end

block_err_r = zeros(length(ns),length(epsilons));
ber_r = zeros(length(ns),length(epsilons));
numofbits = 1e5;
for i = 1:length(ns)
    fprintf('\nn = %d\n',ns(i))
        Lf = [0;ns(i)];
        Rf = [0;0;0;0;0.4*ns(i)];
        H = poly2mat(Rf,Lf);
        H(:,~sum(H)) = [];
        n = size(H,2);
        blocks = ceil(numofbits/n);
        data = zeros(blocks,n);%encoded data
    for j = 1:length(epsilons)
        tic
        fprintf('e = %d\n',epsilons(j))
        channel_data = data;%simulating bec
        add_er = rand(100,n);
        channel_data(add_er<epsilon) = nan;
        corrected = bp_dec(channel_data,H,epsilons(j));
        %Calculating Performance
        iscor = zeros(size(data,1),1);
        corr_bits = iscor;
        for k = 1:size(data,1)
            corr_bits(k) = sum(corrected(k,:)==data(k,:));
            iscor(k) = corr_bits(k)==n;
        end
        block_err_r(i,j) = 1-sum(iscor)/length(iscor);
        ber_r(i,j) = 1-sum(corr_bits)/numel(data);
        toc
    end
end

plot(max(-6,log10(ber)'))
xticks(1:6)
xticklabels(string(0.05:0.05:0.3))
yticks(-6:0)
ylim([-6,0])
ticking = "10^{" + string(-6:0)+"}";
ticking(1) = 0;
ticking(end) = 1;
yticklabels(ticking)
ylabel('Bit Error Rate')
xlabel('erasure probability (e)')
hold on
colors = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30"];
for i = 1:5
    plot(max(-6,log10(ber_r(i,:)')),'color',colors(i),'LineStyle','--')
end
grid on
legend(["n = "+string([100,1000,5000,10000,20000]+" (irregular)"),"n = "+string([100,1000,5000,10000,20000]+" (regular)")])
title("Performance of LDPC")
