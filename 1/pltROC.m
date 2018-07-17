function plt = pltROC(peakRatio)
thresHold = 1.4:1e-4:1.6;
trueSat = 14;
numMonteCarlo = 1000;
numSat = 37;
p_fa = zeros(length(peakRatio),length(thresHold));
p_d  = p_fa;
for ii = 1:length(peakRatio)
    for jj = 1:length(thresHold)
        [a,~] = ind2sub([numSat,numMonteCarlo],find(peakRatio{ii}>thresHold(jj)));
        trueDete  = sum(a == trueSat);
        falseDete = length(a) - trueDete;
        p_fa(ii,jj) = falseDete/numMonteCarlo/36;
        p_d(ii,jj)  = trueDete/  numMonteCarlo;
    end
    idx = find(p_fa(ii,:)<1e-3);
    plt(ii) = p_d(ii,idx(1));
end
end
