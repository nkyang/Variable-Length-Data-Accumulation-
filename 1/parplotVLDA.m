fileIdx{1} = 36:61;
fileIdx{2} = 36:51;
fileIdx{3} = 31:51;
fileIdx{4} = 26:46;
fileIdx{5} = 26:41;
fileIdx{6} = 23:38;
fileIdx{7} = 21:36;
fileIdx{8} = 21:36;
fileIdx{9} = 18:33;
fileIdx{10} = 13:26;
fileList = dir('F:\BDSsim\1\simData\BDSsim_4M_*dB.bin');
k1 = 0;
k2 = 0;
for ii = 1:numel(fileIdx)
    k2 = k1 + length(fileIdx{ii});
    ROC = pltROC(peakRatio(k1+1:k2));
    clear CN0;
    for jj = k1+1:k2
        CN0(jj-k1) = str2double(fileList(fileIdx{ii}(jj-k1)).name((end-8):(end-6)))/10;
    end
    plot(CN0,ROC);
    hold on;
    k1 = k2;
end