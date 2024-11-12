function ecc_plot_1(data,axes,modulation_ind,indi1,indi2)
    M = data.M;
    SNR = data.SNR;
    d = data.d;
    data_rate = data.data_rate;
    mod_name = data.mod_name;
    data = data.err;
    
    if strcmp(axes{1},'bit')
        y = 1;
        yname = 'Bit Error Rate';yticking = cell(10,1);
        yticking{1} = '0';
        yticking{10} = '1';
        zer = -9;
        for i= 1:8
            yticking{10-i} = ['10^{',num2str(-i),'}'];
        end
    elseif strcmp(axes{1},'block')
        y = 2;
        yname = 'Block Error Rate';
        yticking = cell(8,1);
        yticks(-7:0)
        ylim([-7,0])
        yticking{1} = '0';
        yticking{8} = '1';
        zer = -7; 
        for i= 1:6
            yticking{8-i} = ['10^{',num2str(-i),'}'];
        end
    else
        error('Wrong metric name.')
    end
    

    if strcmp(axes{2},'d')
        xname = 'codeword length (n)';
        xticking = string(2.^d-1);
        if strcmp(axes{3},'SNR')
            zname = "SNR = "+string(10*log10(SNR(indi1)))+" dB";
            f = figure('WindowState','maximized');
            plot(max(zer,log10(reshape(data(indi2,:,indi1(1),modulation_ind,y),length(d),1))))
            hold on
            for i = 2:length(indi1)
                plot(max(zer,log10(reshape(data(indi2,:,indi1(i),modulation_ind,y),length(d),1))))
            end
            legend(zname)
            xlabel(xname)
            ylabel(yname)
            yticklabels(yticking)
            xticks(d-1)
            xticklabels(xticking)
            tname = ['M = ',num2str(M(indi2)),' (',mod_name{modulation_ind},')'];
            title(tname)
        elseif strcmp(axes{3},'M')
            zname = 'M = '+ string(M(indi1));
            f = figure('WindowState','maximized');
            plot(max(zer,log10(reshape(data(indi1(1),:,indi2,modulation_ind,y),length(d),1))))
            hold on
            for i = 2:length(indi1)
                plot(max(zer,log10(reshape(data(indi1(i),:,indi2,modulation_ind,y),length(d),1))))
            end
            legend(zname)
            xlabel(xname)
            ylabel(yname)
            yticklabels(yticking)
            xticks(d-1)
            xticklabels(xticking)
            tname = ['SNR = ',num2str(10*log10(SNR(indi2))),' dB (',mod_name{modulation_ind},')'];
            title(tname)
        else
            error('Wrong variable name.')
        end
    elseif strcmp(axes{2},'SNR')
        xname = 'SNR';
        xticking = string(10*log10(SNR));
        if strcmp(axes{3},'d')
            zname = "n = "+string(2.^d(indi1)'-1);
            f = figure('WindowState','maximized');
            plot(max(zer,log10(reshape(data(indi2,indi1(1),:,modulation_ind,y),length(SNR),1))))
            hold on
            for i = 2:length(indi1)
                plot(max(zer,log10(reshape(data(indi2,indi1(i),:,modulation_ind,y),length(SNR),1))))
            end
            legend(zname)
            xlabel(xname)
            ylabel(yname)
            yticklabels(yticking)
            xticklabels(xticking)
            tname = ['M = ',num2str(M(indi2)),' (',mod_name{modulation_ind},')'];
            title(tname)
        elseif strcmp(axes{3},'M')
            zname = 'M = '+ string(M(indi1));
            f = figure('WindowState','maximized');
            plot(max(zer,log10(reshape(data(indi1(1),indi2,:,modulation_ind,y),length(SNR),1))))
            hold on
            for i = 2:length(indi1)
                plot(max(zer,log10(reshape(data(indi1(i),indi2,:,modulation_ind,y),length(SNR),1))))
            end
            legend(zname)
            xlabel(xname)
            ylabel(yname)
            yticklabels(yticking)
            xticklabels(xticking)
            tname = ['n = ',num2str(2.^d(indi2)-1),' (',mod_name{modulation_ind},')'];
            title(tname)
        else
            error('Wrong variable name.')
        end
    elseif strcmp(axes{2},'M')
        x = 1;
        xname = 'modulation order (M)';
        xticking = string(M);
        if strcmp(axes{3},'d')
            zname = "n = "+string(2.^d(indi1)'-1);
            f = figure('WindowState','maximized');
            plot(max(zer,log10(reshape(data(:,indi1(1),indi2,modulation_ind,y),length(M),1))))
            hold on
            for i = 2:length(indi1)
                plot(max(zer,log10(reshape(data(:,indi1(i),indi2,modulation_ind,y),length(M),1))))
            end
            legend(zname)
            xlabel(xname)
            ylabel(yname)
            yticklabels(yticking)
            xticks(1:length(M))
            xticklabels(xticking)
            tname = ['SNR = ',num2str(10*log10(SNR(indi2))),' dB (',mod_name{modulation_ind},')'];
            title(tname)
        elseif strcmp(axes{3},'SNR')
            zname = 'SNR = '+ string(10*log10(SNR(indi1)));
            f = figure('WindowState','maximized');
            plot(max(zer,log10(reshape(data(:,indi2,indi1(1),modulation_ind,y),length(M),1))))
            hold on
            for i = 2:length(indi1)
                plot(max(zer,log10(reshape(data(:,indi2,indi1(i),modulation_ind,y),length(M),1))))
            end
            legend(zname)
            xlabel(xname)
            ylabel(yname)
            yticklabels(yticking)
            xticks(1:length(M))
            xticklabels(xticking)
            tname = ['n = ',num2str(2.^d(indi2)-1),' (',mod_name{modulation_ind},')'];
            title(tname)
        else
            error('Wrong variable name.')
        end
    else
        error('Wrong variable name.')
    end
    grid on
    xlim([1,13])
%     saveas(f,[tname,'.png'])
    

