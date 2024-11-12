%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% parameter_sweep (grid search) for the requested inputs

% Returns bit error rate, block error rate and number of bits used for each
% combination of M,d and modul_name parameter in a 5-D array (M*d*m_n*SNR*err) 
% It also returns the data rate given a symbol time Ts (Default is 1) for
% each combination of M and d in a length(M)*length(d) array.
%
%-IN 
% SNR  scalar
% M,d  vectors
% modul_name cell_array (eg.{'PAM','QAM','PSK'})
% Ts scalar (symbol time)
% nob integer (double) (number of bits for each simulation)
%Ts scalar
%-OUT
% err(5D-array) storing the output of hamming_sim() for each combination of
% the input parameters size =
% [length(M),length(d),length(SNR),length(modul_name),3] last dim is bit-er
% block-er , numofbits
% data_rate(scalar) is the data rate for the given M and Ts size =
% [length(d),length(M)]

function [err,data_rate] = parameter_sweep(d,M,modul_name,SNR,Ts,nob)
    if nargin < 5
        Ts = 1;
        nob = 1e8;
    elseif nargin < 6
        nob = 1e8;
    end
    data_rate = (1-d./(2.^d-1))'.*log2(M)./Ts;%bits/sec
    err = zeros(length(M),length(d),length(modul_name),length(SNR),3);
    for m = 1:length(M)
        fprintf('M = %d\n',M(m))
        tic
        for D = 1:length(d) 
            fprintf('d = %d\n',d(D))
            for MN = 1:length(modul_name) 
                for snr = 1:length(SNR)
                    [err(m,D,snr,MN,1),err(m,D,snr,MN,2),err(m,D,snr,MN,3)] = hamming_sim(SNR(snr),M(m),modul_name{MN},d(D),nob,4e6);%numofbits = 1e8 approx
                end
            end
        end
        toc
    end
