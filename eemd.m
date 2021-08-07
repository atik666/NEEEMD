% y=xlsread('some.data');       % load a signal.
% num_IMF = 5;                      % numbers of IMF
% NR = 10;                      % value of ensemble
% Nstd = 0.3;                   % param to white noise
% 
% IMF1=eemd(y,num_IMF,NR,Nstd);

function [modes,residual] = eemd(y, num_IMF, NR, NstdMax, NstdMin)
stdy = std(y);
if stdy < 0.01
    stdy = 1;
end
y = y ./ stdy;
siz = length(y);
modes = zeros(siz,num_IMF);
res = zeros(siz,1);

for k = 1:NR
    disp(['Ensemble number #' num2str(k)]);

    Nstd = (NstdMax-NstdMin).*rand(1,1) + NstdMin; % Generating random std of white noise

    x = randn(siz,1);
    x = x - mean(x); % genrating 0 mean and 1 std
    x = x -std(x);
    
    wn{k} = x.*Nstd;
    
    y1 = y + wn{k};
    [imf, resN] = emd(y1,'MaxNumIMF',num_IMF); % obtaining IMFs and residuals
    modes = modes + imf;
    res = res + resN; 
end
modes = modes .* stdy ./ (NR); % Final IMFs
residual = res ./ NR; % Final residuals
end
