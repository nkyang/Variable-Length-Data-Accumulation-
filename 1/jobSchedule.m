function peakRatio = jobSchedule
c = parcluster('JB1');
job = createJob(c,'AttachedFiles',{'F:\BDSsim\1\vlda_par.m','F:\BDSsim\1\vlda.mexw64'});
% work_num = 12;
fileList = dir('F:\BDSsim\1\simData\BDSsim_4M_*dB.bin');

% load F:\BDSsim\1\MonterCarlo.mat rStart
% rStart = randi(4.9e7,1,1000);
% save F:\BDSsim\1\newStart.mat rStart
load F:\BDSsim\1\newStart.mat rStart
% interTime = 1e3*[1,2,3,4,5,6,7,8,10,16];
% fileIdx{1} = 36:61;
% fileIdx{2} = 36:51;
% fileIdx{3} = 31:51;
% fileIdx{4} = 26:46;
% fileIdx{5} = 26:41;
% fileIdx{6} = 23:38;
% fileIdx{7} = 21:36;
% fileIdx{8} = 21:36;
% fileIdx{9} = 18:33;
% fileIdx{10} = 13:26;
interTime = 1e3*8;
fileIdx{1} = 21:36;

sub_list = cell(numel(cell2mat(fileIdx)),3);  
kk = 0;
for ii = 1:length(interTime)
    for jj = 1:length(fileIdx{ii})
        kk = kk + 1;
        sub_list{kk,1} = interTime(ii);
        sub_list{kk,2} = fileIdx{ii}(jj);
        sub_list{kk,3} = rStart;
        createTask(job, @vlda_par, 1, {sub_list(kk,:)});
    end
end
submit(job);
wait(job);
peakRatio = fetchOutputs(job);
% plt = pltROC(peakRatio);
end

    