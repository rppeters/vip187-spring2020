[xx,fs] = audioread('a.wav');
tt=(0:length(xx)-1)/fs; %tt is time in s 

middle = length(xx)*fs/2;
duration = length(xx)*fs/2;
fsk = fs;

%apply lowpass under 5kHz
Fn = fs/2;                                                      % Nyquist Frequency (Hz)
Wp = 0.0100;                                                     % Passband Frequency For Lowpass Filter (Hz)
Ws = 0.5000;                                                    % Stopband Frequency For Lowpass Filter (Hz)
Rp =  1;                                                        % Passband Ripple For Lowpass Filter (dB)
Rs = 50;                                                        % Stopband Ripple (Attenuation) For Lowpass Filter (dB)
[n,Wp] = ellipord(Wp,Ws,Rp,Rs);                                 % Calculate Filter Order
[z,p,k] = ellip(n,Rp,Rs,Wp);                                    % Calculate Filter
[sos,g] = zp2sos(z,p,k);                                        % Second-Order-Section For Stability
%xx = filtfilt(sos,g,xx);                                     % Filter Signal

%full spectrogram
[B,F,T,M,xx1,tt1] = myspec(middle, duration, xx, tt, fsk,5000,1);

%ei vowel
[B2, F2, T2, M2, xx2, tt2] =myspec(0.47,0.13,xx,tt,fsk, 5000, 2);
xs = length(xx2);

b= floor(xs/5); %beginning
m= floor(xs/2); %middle
e= floor(xs*3/4); %end
wsize = 4000;

%apply window to time series
w1=[0*(1:b-wsize/2)'; hann(wsize); 0*(1:xs-b-(wsize/2))']; %beginning vowel
length(xx2)
length(w1)
z1=xx2.*w1;
w2=[0*(1:m-wsize/2)'; hann(wsize); 0*(1:xs-m-(wsize/2))']; %middle vowel
z2=xx2.*w2;
w3=[0*(1:e-wsize/2)'; hann(wsize); 0*(1:xs-e-(wsize/2))']; %end vowel
z3=xx2.*w3;

% get FFT spectra
y1=abs(fft(z1));
y2=abs(fft(z2));
y3=abs(fft(z3));
% construct frequency vector
fr=(0:length(y1)-1)*fsk/length(y1);
ifr=fr<fsk/2;

% obtain frequency of spectral maxima
frmax1=fr(find(y1==max(y1),1,'first'));
frmax2=fr(find(y2==max(y2),1,'first'));
frmax3=fr(find(y3==max(y3),1,'first'));

%plot variables
xbounds = [1 4000];

%plot time series and spectra side by side
figure(2)
set(gcf,'position',[100, 100, 560, 600]) 
subplot(321)
hp=plot(tt2,xx2,'k',tt2,z1,'g'); xlim(tt2([1 end]))
set(hp(2),'linewidth',2)
ylabel('rel. Amplitude')
subplot(322)
plot(fr(ifr),y1(ifr),'k'); xlim(xbounds)
text(100,max(y1),sprintf('max = %.1f kHz',frmax1))
subplot(323)
hp=plot(tt2,xx2,'k',tt2,z2,'g'); xlim(tt2([1 end]))
set(hp(2),'linewidth',2)
ylabel('rel. Amplitude')
subplot(324)
plot(fr(ifr),y2(ifr),'k'); xlim(xbounds)
text(100,max(y2),sprintf('max = %.1f kHz',frmax2))
subplot(325)
hp=plot(tt2,xx2,'k',tt2,z3,'g'); xlim(tt2([1 end]))
set(hp(2),'linewidth',2)
xlabel('Time [ms]')
ylabel('rel. Amplitude')
subplot(326)
plot(fr(ifr),y3(ifr),'k');xlim(xbounds)
text(100,max(y3),sprintf('max = %.1f kHz',frmax3))
xlabel('Frequency [kHz]')

%overlay spectra
figure(3)
xlabel('Time [ms]')
ylabel('rel. Amplitude')
plot(fr(ifr),y1(ifr),'b');
hold on
plot(fr(ifr),y2(ifr),'g');
plot(fr(ifr),y3(ifr),'r');xlim(xbounds)
hold off
text(100,max(y3),sprintf('max = %.1f kHz',frmax3))
xlabel('Frequency [kHz]')
legend('Beginning', 'Middle', 'End')
title('Freq Change in [ei] overtime')

%get spectrogram for specific start time and duration (s,d)
function [B,F,T,M,xx1,tt1] = myspec(start, dur, xx, tt, fsk, nfft, fignum)
    %limit to region
    ts=start+[-1 1]*dur; 
    tsel=tt>ts(1) & tt<ts(2); %time selected
    tt1=(tt(tsel)-start)*1000;
    xx1=xx(tsel); %reduce time-series to just selected range

    %reduce vertical-offset???
    xx1=xx1-mean(xx1);

    [B,F,T] = specgram(xx1,nfft,fsk,hann(32),12); 
    M=20*log10(abs(B));

    %plot spectrogram
    figure(fignum)
    imagesc(T-dur,F,M); axis xy
    grid on
    colormap(1-gray(12))
    cl=caxis;
    caxis(-20+[-60 0])
    %caxis(cl(2)+[-60 0])
    hb=colorbar;
    set(get(hb,'title'),'string','[dB]')
    
    ylim([0 6*1000])
    xlabel('Time [ms]')
    ylabel('Frequency [Hz]')
    title('Bake [beik]')
end


