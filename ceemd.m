% y=xlsread('some.data');       % load a signal.
% aim = 5;                      % numbers of IMF
% NR = 10;                      % value of ensemble
% Nstd = 0.3;                   % param to white noise
% 
% IMF1=eemd(y,aim,NR,Nstd);

function [modes,residual] = ceemd(y, num_IMF, NR, NstdMax, NstdMin)
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
    y2 = y - wn{k};
    
    [imf1, res1] = emd(y1,'MaxNumIMF',num_IMF); % obtaining IMFs and residuals
    modes = modes + imf1;
    res = res + res1; 
    if Nstd > 0 && NR > 1
        [imf2, res2] = emd(y1,'MaxNumIMF',num_IMF); % obtaining IMFs and residuals
        modes = modes + imf2;
        res = res + res2; 
    end

end
modes = modes .* stdy ./ (NR);
if Nstd > 0 && NR > 1
    modes = modes ./ 2;
end
residual = (res ./ NR)./ 2 ; % Final residuals
end