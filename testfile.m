%Task 1
%This script generates binary data, encodes them, simulates a binary
%symetric channel and then decodes them and calculates block error rate and
%ber.
%
%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418
close all
clear
clc
%% Init
rng('default')
SNR = 10;
M =2;
[H,G,n,k] = hammgen(12);
modul_name = 'PAM';
p = pb_err(M,modul_name,SNR);
%p = 0.1;
% k = 4;
% n = 7;
%calculating batch length so that the script doesnt utilize more than 
% 4GB of RAM 
l = floor(4e9/(4*(2*n-k)*M*k*8));%we need 2 copies of the data and some calculations
%add a copy temporarily (pretty good estimation after some simulations)
%l = 1e4;%batch size
if l < 1
    error('Not enough memory.')
end
bnum = 1;%batches
data_sz = l*M*n*k; 
blockerr_c = 0;%block errors counter
biterr_c = 0;%bit errors counter
e = sparse(eye(n));%Hamming only
s = e*H';
for j = 1:bnum
    batchtime = tic;
    fprintf('Batch num %d\nBatch size %d blocks:\n\n',j,l*M*n)
    disp('Data Generation..')
    tic
    data = randi([0,1],data_sz,1,'double');%randi([0,255],data_sz,1,'uint8');
    data = reshape(data,k,[]);
    toc
    disp('Encoding...')
    tic
    temp = mod(data'*G,2)';
    toc
    disp('Transmitting...')
    tic
    temp = bsc(temp,p);
    toc
    disp('Syndrome...')
    tic
    syn = mod((H*temp)',2);
    toc
    disp('Error correcting...')
    tic
    for i =1:2^(n-k)-1
        ind = sum(syn == s(i,:),2);
        temp(i,ind==n-k) = mod(temp(i,ind==n-k) + 1,2);%hamming only
    end
    toc
    disp('Decoding...')
    tic
    temp = temp(n-k+1:end,:);
    toc
    disp('Error calculations...')
    tic
    %bter = sum(~(temp==data));
    biterr_c = sum(sum(~(temp==data)));
    blockerr_c = sum(logical(sum(~(temp==data))));
    toc
    disp('Batch total time.')
    toc(batchtime)
end
biterr_c = biterr_c/(data_sz);
blockerr_c = blockerr_c/(data_sz);





