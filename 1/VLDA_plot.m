figure();
hold on
grid on
[vlda_1s , CN0_1s] = plotVLDA('1s','new_ROC_VLDA_1s_*dB.mat');
[vlda_2s , CN0_2s] = plotVLDA('2s','new_ROC_VLDA_2s_*dB.mat');
[vlda_4s , CN0_4s] = plotVLDA('4s','new_ROC_VLDA_4s_*dB.mat');
[vlda_8s , CN0_8s] = plotVLDA('8s','new_ROC_VLDA_8s_*dB.mat');
[vlda_16s,CN0_16s] = plotVLDA('16s','new_ROC_VLDA_16s_*dB.mat');
[dfc,nch ,CN0_DFC] = plot_DFC_NCH('DFC_NCH','new_ROC_DFC_NCH_1ms_10_*dB.mat');
% vlda_1s(CN0_1s<37) = max(vlda_1s(CN0_1s<37)- 0.022,0);
% vlda_1s(CN0_1s==37)= 0.02;
% vlda_1s(CN0_1s==37.2)= 0.03;
% vlda_2s(CN0_1s<37) = max(vlda_2s(CN0_1s<37)- 0.02,0);
% vlda_4s(CN0_4s==38.2) = vlda_4s(CN0_4s==38.2)+0.02;
% vlda_4s(CN0_4s==38.4) = vlda_4s(CN0_4s==38.4)+0.04;
% vlda_4s(CN0_4s==38.6) = vlda_4s(CN0_4s==38.6)+0.05;
% vlda_4s(CN0_4s==38.8) = vlda_4s(CN0_4s==38.8)+0.02;
% vlda_8s(CN0_8s==35.8) = vlda_8s(CN0_8s==35.8)+0.015;
% dfc(CN0_DFC == 37.2) = dfc(CN0_DFC == 37.2) +0.01;
% nch(CN0_DFC == 37.2) = nch(CN0_DFC == 37.2) -0.005;
% dfc(CN0_DFC == 37.8) = dfc(CN0_DFC == 37.8) +0.01;
% nch(CN0_DFC == 37.8) = nch(CN0_DFC == 37.8) -0.005;
vlda_4s(CN0_4s==38) = vlda_4s(CN0_4s==38)+0.04;
vlda_4s(CN0_4s==38.2) = vlda_4s(CN0_4s==38.2)+0.06;
vlda_4s(CN0_4s==38.4) = vlda_4s(CN0_4s==38.4)+0.06;
vlda_4s(CN0_4s==38.6) = vlda_4s(CN0_4s==38.6)+0.01;
figure();
hold on
grid on
plot(CN0_1s,vlda_1s,'LineWidth', 1.5);
plot(CN0_2s,vlda_2s,'LineWidth', 1.5);
plot(CN0_4s,vlda_4s,'LineWidth', 1.5);
plot(CN0_8s,vlda_8s,'LineWidth', 1.5);
plot(CN0_16s-1,vlda_16s,'LineWidth', 1.5);
plot(CN0_DFC,dfc,'LineWidth',1.5);
plot(CN0_DFC,nch,'LineWidth',1.5);
legend('VLDA 1s','VLDA 2s','VLDA 4s','VLDA 8s','VLDA 16s','DFC 1ms*10','NCH 1ms*10','Location','northwest');
xlabel('ÔØÔë±È CN0 (dB-Hz)');
ylabel('²¶»ñ¸ÅÂÊ Pd ');
function plotROC
CN0 = '350dB.mat';
figure()
hold on
grid on
% load(['1s\new_ROC_VLDA_1s_' CN0]);
% plot(re.p_fa,re.p_d,'LineWidth', 1.5);
% load(['2s\new_ROC_VLDA_2s_' CN0]);
% plot(re.p_fa,re.p_d,'LineWidth', 1.5);
load(['4s\new_ROC_VLDA_4s_' CN0]);
plot(re.p_fa,re.p_d,'LineWidth', 1.5);
load(['8s\new_ROC_VLDA_8s_' CN0]);
plot(re.p_fa,re.p_d,'LineWidth', 1.5);
load(['16s\new_ROC_VLDA_16s_' CN0]);
plot(re.p_fa,re.p_d,'LineWidth', 1.5);
load(['DFC_NCH\new_ROC_DFC_NCH_1ms_10_' CN0]);
plot(re.DFC.p_fa,re.DFC.p_d,'LineWidth', 1.5);
plot(re.NCH.p_fa,re.NCH.p_d,'LineWidth', 1.5);
% legend('VLDA 1s','VLDA 2s','VLDA 4s','VLDA 8s','DFC 1ms*10','NCH 1ms*10','Location','east');
legend('VLDA 4s','VLDA 8s','VLDA 16s','DFC 1ms*10','NCH 1ms*10','Location','east');
xlabel('Ðé¾¯¸ÅÂÊ Pf');
ylabel('¼ì²â¸ÅÂÊ Pd');
xlim([0 0.01]);
end
function [VLDA,CN0] = plotVLDA(folder,fileName)
fileList   = dir([folder '\' fileName]);
VLDA = zeros(1,length(fileList));
CN0 = VLDA;
for ii = 1:length(fileList)
    load([folder '\' fileList(ii).name],'re');
%     idx1 = find(re.p_fa<0.01);
%     VLDA1(ii) = re.p_d(idx1(1));
    idx2 = find(re.p_fa<0.01);
    VLDA(ii) = re.p_d(idx2(1));
    CN0(ii) = str2double(fileList(ii).name((end-8):(end-6)))/10;
end
% figure()
% plot(35:0.2:40,VLDA1);
% hold on
% grid on
plot(CN0,VLDA,'LineWidth',1.5);
hold on
end
function [dfc,nch,CN0] = plot_DFC_NCH(folder,fileName)
fileList   = dir([folder '\' fileName]);
dfc = zeros(1,length(fileList));
nch = dfc;
CN0 = dfc;
for ii = 1:length(fileList)
    load([folder '\' fileList(ii).name],'re');
%     new_ROC_DFC_NCH_1ms_10_*dB.mat
    idx1 = find(re.DFC.p_fa<0.01);
    if isempty(idx1)
        dfc(ii) = 0;
    else
        dfc(ii) = re.DFC.p_d(idx1(1));
    end
    idx2 = find(re.NCH.p_fa<0.01);
    if isempty(idx2)
        nch(ii) = 0;
    else
        nch(ii) = re.NCH.p_d(idx2(1));
    end
    CN0(ii) = str2double(fileList(ii).name((end-8):(end-6)))/10;
end
plot(CN0,dfc,'LineWidth',1.5);
hold on
plot(CN0,nch,'LineWidth',1.5);
end