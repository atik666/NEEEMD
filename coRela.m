function corr = coRela(Signal, Imf)

stdSig = std(Signal);
stdImf = std(Imf);
tStd = stdSig.*stdImf;

cMat = cov(Signal, Imf);

corr = cMat ./ tStd;

end