clc;clear;

y = csvread('simulatedSignal.csv');

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
tk = tesk(y);

TESK = [tk1; tk2; tk3; tk4];
sumTesk = vertcat(tk, sum(TESK,2));

rrmse1 = rms(y- (sum(emdIMF,2)+residual1))./rms(sum(emdIMF,2)+residual1);
rrmse2 = rms(y- (sum(eemdIMF,2)+residual2))./rms(sum(eemdIMF,2)+residual2);
rrmse3 = rms(y- (sum(ceemdIMF,2)+residual3))./rms(sum(ceemdIMF,2)+residual3);
rrmse4 = rms(y- (sum(neeemdIMF,2)+residual4))./rms(sum(neeemdIMF,2)+residual4);

rrmse = vertcat(rrmse1,rrmse2,rrmse3,rrmse4);
%csvwrite('tk.csv',sumTesk);

rmse1 = sqrt(mean((y-sum(emdIMF,2)+residual1).^2));
rmse2 = sqrt(mean((y-sum(eemdIMF,2)+residual2).^2));
rmse3 = sqrt(mean((y-sum(ceemdIMF,2)+residual3).^2));
rmse4 = sqrt(mean((y-sum(neeemdIMF,2)+residual4).^2));

for i = 1:4
    subplot(5,1,i)
    plot(eemdIMF(:,i)) 
    set(gca,'xtick',[])
end
subplot(5,1,5)
plot(residual2(:,:))

