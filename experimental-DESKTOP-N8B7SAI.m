clc;clear;close all;

% data = csvread ('D:\OneDrive - Universiti Malaysia Pahang\Atik_Home\Data Files\fyp-velocity signals\Twist_3.csv');
% y = data(1:2500,:);

data = load('D:\OneDrive - ump.edu.my\Atik_Home\Data Files\Blade Data\data_20150516_20hz\t1_20hz_16.mat');
y = data.Channel_003(1:2500,:);

num_IMF = 4;
NR = 100;
NstdMax = 0.2;
NstdMin = 0.1;

[emdIMF,residual1] = emd(y,'MaxNumIMF',num_IMF);
[eemdIMF,residual2] = eemd(y, num_IMF, NR, NstdMax, NstdMin);
[ceemdIMF,residual3] = ceemd(y, num_IMF, NR, NstdMax, NstdMin);
[neeemdIMF,residual4] = neeemd(y, num_IMF, NR, NstdMax, NstdMin);

[tk1,kc1,kes1] = tesk(emdIMF);
[tk2,kc2,kes2] = tesk(eemdIMF);
[tk3,kc3,kes3] = tesk(ceemdIMF);
[tk4,kc4,kes4] = tesk(neeemdIMF);

cor1 = corR(y,emdIMF);
cor2 = corR(y,eemdIMF);
cor3 = corR(y,ceemdIMF);
cor4 = corR(y,neeemdIMF);

cor = [(cor2); (cor3); (cor4)]';

L = 20;
%[y_final, f_final, kurtIter] = med2d(neeemdIMF(:,1),30,[],[],'valid',1);
[~, f, y_final] = momeda(neeemdIMF(:,1),L,[1 1 1],5:0.1:300,0);
%y_final = vertcat((0*ones(1,L))', y_final);
% [~, ~, emdIMF] = momeda(emdIMF(:,1),500,[1 1 1 1 1],5:0.1:300,1);
% [~, ~, eemdIMF] = momeda(eemdIMF(:,1),500,[1 1 1 1 1],5:0.1:300,1);
% [~, ~, ceemdIMF] = momeda(ceemdIMF(:,1),500,[1 1 1 1 1],5:0.1:300,1);

fs = 12000;
xmin = 0;
xmax = 1000;
% subplot(4,1,1)
% [es1,F1] = envspectrum(y,fs);
% plot(F1,es1.*1./max(es1))
% xlim([xmin xmax])
subplot(3,1,1)
[es2,F2] = envspectrum(eemdIMF(:,1),fs);
plot(F2,es2.*1./max(es2))
xlim([xmin xmax])
subplot(3,1,2)
[es3,F3] = envspectrum(ceemdIMF(:,1),fs);
plot(F3,es3.*1./max(es3))
xlim([xmin xmax])
subplot(3,1,3)
[es4,F4] = envspectrum(y_final,fs);
plot(F4,es4.*1./max(es4))
xlim([xmin xmax])

figure(2)
subplot(3,1,1)
plot(eemdIMF(:,1))
subplot(3,1,2)
plot(ceemdIMF(:,1))
subplot(3,1,3)
plot(y_final)

sp1 = sens(emdIMF);
sp2 = sens(eemdIMF);
sp3 = sens(ceemdIMF);
sp4 = sens(neeemdIMF);

sp = [(sp2); (sp3); (sp4)]';