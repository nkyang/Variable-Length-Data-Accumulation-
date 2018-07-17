function peakRatio = vlda_par(interTime,fileIdx,rStart)
%% Initialization ==================================================
samplingFreq = 10e6;
% interTime = 6e3;
% Find number of samples per spreading code
samplesPerCode = 10000;
samplesPerCodeChip = round(samplesPerCode / 2046);
dopplerCode = (-7e3:500/(interTime/1e3):7e3)/763;
flag = int32(sign(dopplerCode));
interSample= interTime*samplesPerCode;
segmentLen = int32(abs(2*round(samplesPerCode*(102.3)./dopplerCode)-1));
segmentNum = int32(floor(int32(interSample)./segmentLen));
totalLen   = interSample - flag.*segmentNum;
% Frequency Domin of ranging code
fftCode    = longCode(samplingFreq,samplesPerCode,1);
ftc        = reshape(fftCode.',samplesPerCode,1,37);
numMonteCarlo = length(rStart);
fileList   = dir('simData\BDSsim_4M_*dB.bin');
%--- Initialize arrays to speed up the code -------------------------------
peakRatio = ones(37,numMonteCarlo);
%% 
fid  = fopen(fileList(fileIdx).name,'rb');
data = fread(fid,5e7+interSample,'int8=>int16');
fclose(fid);
data = data.*circshift(data,[22 0]);
for jj = 1:numMonteCarlo
    cumSignal = vlda(data(rStart(jj):rStart(jj)+interSample+1*samplesPerCode),segmentLen,totalLen,flag);
    signalFreqDom = fft(cumSignal,samplesPerCode,1);
    results  = abs(ifft(bsxfun(@times,signalFreqDom,ftc),samplesPerCode,1));
    [peakSize,codePhase] = max(max(results,[],2),[],1);
    [~,codeDoppleIndex] = max(max(results,[],1),[],2);
    for kk = 1:37
        id1 = codePhase(kk) - samplesPerCodeChip;
        id2 = codePhase(kk) + samplesPerCodeChip;
        A   = 1:samplesPerCode;
        codePhaseRange = A((A<id1)&(A>id2-samplesPerCode)|((A>id2)&(A<id1+samplesPerCode)));
        peakRatio(kk,jj) = peakSize(kk)/max(results(codePhaseRange,codeDoppleIndex(kk),kk));
    end
end
end

%%
function fftCode = longCode(samplingFreq,samplesPerCode,n)
longCode = makeRangTable(samplingFreq,samplesPerCode,n);
longCodeDelay = longCode.*circshift(longCode,[0,22]);
longCodeDelay = double(longCodeDelay);
fftCode  = conj(fft(longCodeDelay,[],2));
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
