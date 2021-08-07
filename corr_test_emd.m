clear;clc;

base1 = csvread('D:\OneDrive - Universiti Malaysia Pahang\Atik_Home\Data Files\fyp-velocity signals\Baseline_1.csv');

num_IMF = 5;
NR = 100;
NstdMax = 0.2;
NstdMin = 0.1;

data_cycle = 1200;
L = 4;
k=1;
for i=1:L
    E = base1(k:k+floor(data_cycle)-1);
    [nEEEMD_IMFs{i,1}, EEMD_IMFs{i,1}] =  nEEMD(E, num_IMF, NR, NstdMax, NstdMin);
    EMD_IMFs{i,1} = emd(E,'MaxNumIMF',num_IMF);
end  

a1 = {};
a2 = {};
a3 = {};
for i = 1: size(EMD_IMFs{1,1},2)
    a1{1,i} = corrcoef(E,EMD_IMFs{1,1}(:,i));
    a2{1,i} = corrcoef(E,EEMD_IMFs{1,1}(:,i));
    a3{1,i} = corrcoef(E,nEEEMD_IMFs{1,1}(:,i));
end
a1 = cell2mat(a1);
a2 = cell2mat(a2);
a3 = cell2mat(a3);

aaa = vertcat(a1,a2,a3);

for i = 1:2
    tkEMD{i,1} = tesk(EMD_IMFs{i,1});
    tkEEMD{i,1} = tesk(EEMD_IMFs{i,1});
    tknEEEMD{i,1} = tesk(nEEEMD_IMFs{i,1});
    energy{i,1} = kurtosis(nEEEMD_IMFs{i,1});
    energyEEMD{i,1} = kurtosis(EEMD_IMFs{i,1});
    energyEMD{i,1} = kurtosis(EMD_IMFs{i,1});
end

TK = [sum(tkEMD{2,1}), sum(tkEEMD{2,1}), sum(tknEEEMD{2,1})];

coEMD = coRela(base1(1:1200)', EMD_IMFs{1,1}(:,1));
 
aa = (base1(1:1200,:));
 
aa1 = sum(EEMD_IMFs{4,1}, 2);

rms2 = rms(abs(aa-aa1))./rms(aa1);

noise1 = snr(base1(1:1200,1),sum(EMD_IMFs{1,1},2));
noise2 = snr(base1(1:1200,1),sum(EEMD_IMFs{1,1},2));
noise3 = snr(base1(1:1200,1),sum(nEEEMD_IMFs{1,1},2));

power = AV_PSD(base1(1:1200,1));
noise = power - AV_PSD(sum(nEEEMD_IMFs{3,1},2));
%noise1 = 10.*log(power./noise);
noise3 = snr(power,noise);