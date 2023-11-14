%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%IN
% M is the order of the modulation must be positive integer
% modul_name the name of the modulation accepted values are
% ('PAM','PSK','QAM')
% SNR is the SNR of the link (double)
%OUT
%p is the error probability of the bsc (double)

function p = pb_err(M,modul_name,SNR)
    switch modul_name
        case 'PAM'
            p = qfunc((2*SNR)^(.5));
        case 'PSK'
            p = M*qfunc((2*SNR)^(.5))/4;
        case 'QAM'
            p = M*qfunc((2*SNR)^(.5))/2;
        otherwise
            error('You entered unsupported modulation name!')
    end