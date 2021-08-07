function [Power_full_signal, xfft] = AV_PSD(Data)
fs=1000;
nfft=length(Data);
nfft2=2^nextpow2(nfft);
ff=fft(Data,nfft2);
fff=ff(1:nfft2/2,:);
xfft=fs*(0:nfft2/2-1)/nfft2;
Power_full_signal=(1/nfft2)*(abs(fff).^2);
%Av_Power = mean(Power_full_signal); %max=81.1%
end


