clc;clear;

y = csvread ('simulatedSignal.csv');

num_IMF = 4;
NR = 100;
NstdMax = 0.2;
NstdMin = 0.1;

[emdIMF,residual1] = emd(y);
[eemdIMF,residual2] = eemd(y, num_IMF, NR, NstdMax, NstdMin);
[ceemdIMF,residual3] = ceemd(y, num_IMF, NR, NstdMax, NstdMin);
[neeemdIMF,residual4] = neeemd(y, num_IMF, NR, NstdMax, NstdMin);

noise1 = mean(AV_PSD(emdIMF)) ./ mean(AV_PSD(residual1));
noise2 = mean(AV_PSD(eemdIMF)) ./ mean(AV_PSD(residual2));
noise3 = mean(AV_PSD(ceemdIMF)) ./ mean(AV_PSD(residual3));
noise4 = mean(AV_PSD(neeemdIMF)) ./ mean(AV_PSD(residual4));

[tk1,kc1, kes1] = tesk(emdIMF);
[tk2,kc2, kes2] = tesk(eemdIMF);
[tk3,kc3, kes3] = tesk(ceemdIMF);
[tk4,kc4, kes4] = tesk(neeemdIMF);
[tk,kc, kes] = tesk(y);

TESK = [tk1; tk2; tk3; tk4];
sumTesk = sum(TESK,2);

cor1 = corR(y,emdIMF);
cor2 = corR(y,eemdIMF);
cor3 = corR(y,ceemdIMF);
cor4 = corR(y,neeemdIMF);

cor = [sum(cor1); sum(cor2); sum(cor3); sum(cor4)];

fs = 1000;
subplot(4,1,1)
es1 = envspectrum(emdIMF(:,1),fs);
plot(es1)
subplot(4,1,2)
es2 = envspectrum(eemdIMF(:,1),fs);
plot(es2)
subplot(4,1,3)
es3 = envspectrum(ceemdIMF(:,1),fs);
plot(es3)
subplot(4,1,4)
es4 = envspectrum(neeemdIMF(:,2),fs);
plot(es4)

