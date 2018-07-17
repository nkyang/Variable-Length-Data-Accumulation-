function [DFC,NCH] = DFC_NCH_MonteCarlo
%% Initialization ==================================================
samplingFreq = 10e6;
interFreq = 2.5e6;

% Find number of samples per spreading code
samplesPerCode = 10000;
samplesPerCodeChip = round(samplesPerCode / 2046);
% Find sampling period
ts = 1 / samplingFreq;
% Find phase points of the local carrier wave
phasePoints_D2 = (0 : (11*samplesPerCode-1)) * 2 * pi * ts;
% Number of the frequency bins for the given acquisition band (500Hz steps)
% Doppler shift +/- 7KHz
numberOfFreqBins = round(14 * 2) + 1; 
% Generate all ranging code and sample them according to the sampling freq.
rangCodesTable  = makeRangTable(samplingFreq,samplesPerCode,1)';
% Frequency Domin of ranging code
rangCodeFreqDom = reshape(conj(fft(rangCodesTable,samplesPerCode,1)),[samplesPerCode,1,1,37]);
% Carrier frequencies of the frequency bins
frqBins   = interFreq - (14/2) * 1000 + 500 * (0:numberOfFreqBins - 1);
carr = exp(1i*bsxfun(@times,frqBins',phasePoints_D2).');
carr = reshape(carr,samplesPerCode,11,numberOfFreqBins);
% trueSat    = [1,2,3,4,5,6,7,9,13,14,17,32];
trueSat = 14;
wrongSat   = setdiff(1:37,trueSat);
numMonteCarlo = 1000;
% thresHold  = 1.2:0.001:1.7;
% thresHold  = 1.2:0.0002:1.35;
thresHold = 1680:0.5:1920;
sqThresHold = thresHold.^2;
fileList  = dir('BDSsim_4M_*dB.bin');
%--- Initialize arrays to speed up the code -------------------------------
DFC.nFalse = zeros(1,length(thresHold));
DFC.nCorre = zeros(length(trueSat),length(thresHold));
DFC = repmat(DFC,[1,length(fileList)]);
NCH = DFC;
A   = 1:samplesPerCode;
load MonterCarlo.mat rStart
for ii = 53:54
    fid = fopen(fileList(ii).name,'rb');
    longSignal = fread(fid,5e7,'int8=>int8');
    fclose(fid);
    for jj = 1:numMonteCarlo
        signal   = double(reshape(longSignal(rStart(jj):11*samplesPerCode+rStart(jj)-1),samplesPerCode,11));
        acqRes   = ifft(bsxfun(@times,fft(carr.*signal,samplesPerCode,1),rangCodeFreqDom),samplesPerCode,1);
        %% DFC
        results  = sum(abs(acqRes(:,1:10,:,:).*acqRes(:,2:11,:,:)),2);
%         [peakSize,codeDoppleIndex] = max(max(results,[],1),[],3);
%         [~,codePhase] =max(max(results,[],3),[],1);
%         for kk = 1:37
%             id1 = codePhase(kk) - samplesPerCodeChip;
%             id2 = codePhase(kk) + samplesPerCodeChip;
%             codePhaseRange = A((A<id1)&(A>id2-samplesPerCode)|((A>id2)&(A<id1+samplesPerCode)));
%             peakRatio(kk,1) = peakSize(kk)/max(results(codePhaseRange,1,codeDoppleIndex(kk),kk));
%         end
%         decision = bsxfun(@ge,peakRatio,sqThresHold);
        peakSize = realsqrt(squeeze(max(max(results,[],1),[],3))/10);
        decision = bsxfun(@ge,peakSize,thresHold);
        DFC(ii).nCorre = DFC(ii).nCorre + decision(trueSat,:);
        DFC(ii).nFalse = DFC(ii).nFalse + sum(decision(wrongSat,:));
        %% NCH
        results = sum(abs(acqRes(:,1:10,:,:)),2);
%         [peakSize,codeDoppleIndex] = max(max(results,[],1),[],3);
%         [~,codePhase] =max(max(results,[],3),[],1);
%         for kk = 1:37
%             id1 = codePhase(kk) - samplesPerCodeChip;
%             id2 = codePhase(kk) + samplesPerCodeChip;
%             codePhaseRange = A((A<id1)&(A>id2-samplesPerCode)|((A>id2)&(A<id1+samplesPerCode)));
%             peakRatio(kk,1) = peakSize(kk)/max(results(codePhaseRange,1,codeDoppleIndex(kk),kk));
%         end
%         decision = bsxfun(@ge,peakRatio,thresHold);
        peakSize = squeeze(max(max(results,[],1),[],3)/10);
        decision = bsxfun(@ge,peakSize,thresHold);
        NCH(ii).nCorre = NCH(ii).nCorre + decision(trueSat,:);
        NCH(ii).nFalse = NCH(ii).nFalse + sum(decision(wrongSat,:));
    end
    DFC(ii).p_fa = DFC(ii).nFalse/(numMonteCarlo*(37 - length(trueSat)));
    DFC(ii).p_d  = DFC(ii).nCorre/ numMonteCarlo;
    NCH(ii).p_fa = NCH(ii).nFalse/(numMonteCarlo*(37 - length(trueSat)));
    NCH(ii).p_d  = NCH(ii).nCorre/ numMonteCarlo;
    re.NCH = NCH(ii);
    re.DFC = DFC(ii);
    eval(['save' ' new_ROC_DFC_NCH_1ms_10_' fileList(ii).name(11:15) '.mat' ' re']);
end
% save ROC_DFC_1ms_10.mat DFC
% save ROC_NCH_1ms_10.mat NCH
end

%% make ranging code Table
function rangCodesTable = makeRangTable(samplingFreq,samplesPerCode,n)
%--- Find number of samples per spreading code ----------------------------
len = round(n * samplesPerCode);
codeFreq = 2.046e6;
%--- Find time constants --------------------------------------------------
ts = 1/samplingFreq;   % Sampling period in sec
tc = 1/codeFreq;  % code chip period in sec
%--- Prepare the output matrix to speed up function -----------------------
rangCodesTable = zeros(37, len);
%% === For all satellite PRN-s ...
for i=1:37
    %--- Generate code for given PRN --------------------------------------
    rangCode = generateRangCode(i);
    rangValueIndex = ceil((ts * (1:len))/tc);
    %--- Correct the last index (due to number rounding issues) -----------
    rangValueIndex = mod(rangValueIndex-1,2046)+1;
    rangValueIndex(end) = 2046;
    %--- Make the digitized version of the ranging code -------------------
    rangCodesTable(i, :) = rangCode(rangValueIndex);   
end % for PRN = 1:37
end


%% Generate ranging code
function rangcode = generateRangCode(PRN)
g2s =  [1335, 466, 633, 497,1466,1276, 736,1004, 498,1688,...
        1337, 468, 499, 944,1468,1278,1689,1338, 636, 500,...
         945,1469,1690, 470, 637, 501, 946,1340, 471, 638,...
         502,1693,1342, 473,1694,1343,1695];
%% Generate g1
g1=zeros(1,2047);
reg=[-1,1,-1,1,-1,1,-1,1,-1,1,-1];
for i=1:2047
    g1(i)       = reg(11);
    saveBit     = -(reg(1)*reg(7)*reg(8)*reg(9)*reg(10)*reg(11));
    reg(2:11)   = reg(1:10);
    reg(1)      = saveBit;
end
%% Generate g2
g2=zeros(1,2047);
reg=[-1,1,-1,1,-1,1,-1,1,-1,1,-1];
for i=1:2047
    g2(i)       = reg(11);
    saveBit     = -(reg(1)*reg(2)*reg(3)*reg(4)*reg(5)*reg(8)*reg(9)*reg(11));
    reg(2:11)   = reg(1:10);
    reg(1)      = saveBit;
end
% g2 offset
g2 = circshift(g2,[0,g2s(PRN)]);
%% Generate ranging code
% 注意需要截短一个码片
rangcode = -(g1(1:2046).*g2(1:2046));
end
