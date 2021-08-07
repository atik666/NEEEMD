function [tk,kc, kes] = tesk(imf)

uc = mean(imf);

expect = (imf - uc).^4;
expect = mean(expect);
mc = std(imf);
mc = mc.^4;
kc = expect./mc;

%esf = abs(fft(hilbert(imf)));
esf = pwelch(imf);

ues = mean(esf);
expes = esf - ues;
expes = (expes).^4;
expes = mean(expes);
mes = std(esf);
mes = (mes).^4;
kes = expes./mes;

tk = kc.*kes;

end