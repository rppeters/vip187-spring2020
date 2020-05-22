[xx,fs] = audioread('a.wav');
tt=(0:length(xx)-1)/fs; %tt is time in s 
fsk = fs;

%reduce to [a] vowel
middle = 0.55;
dur=0.13;
ts=middle+[-1 1]*dur; 
tsel=tt>ts(1) & tt<ts(2); %time selected
tt1=(tt(tsel)-middle)*1000;
xx1=xx(tsel); %reduce time-series to just selected range
xs = length(xx2);

b=floor(xs/2); %middle
wsize = floor(4096*2);

%apply window to time series
w1=[0*(1:b-wsize/2)'; hann(wsize); 0*(1:xs-b-(wsize/2))']; %beginning vowel
z1=xx2.*w1;

% get FFT spectra
y1=abs(fft(z1));

% construct frequency vector
fr=(0:length(y1)-1)*fsk/length(y1);
ifr=fr<fsk/2;

% obtain frequency of spectral maxima
frmax1=fr(find(y1==max(y1),1,'first'));
[pks,locs] = findpeaks(y1,'MinPeakDistance',25);

%plot variables
xbounds = [1 4000];
plot(fr(ifr),y1(ifr),'k'); xlim(xbounds)
hold on
for i=1:100
    plot(locs(i)*3.8,pks(i),'r*');
    text(locs(i)*3.8,pks(i),sprintf('<-- H%d=%f',i,locs(i)*3.8),'FontSize',8,'Rotation',90);
end
hold off

c = 350; %m/s
%speed,lambda,w,T,f,P calculations
f0 = locs(1)*3.8; %fundamental freq from [a]
w = f0*2*pi;
T = 1 / f0;
lambda = c/f0;
k = 2*pi/lambda;
fprintf('Sound Speed=%d m/s \nFreq=%f2 Hz \nW=%f4  \nPeriod=%f4 s\nWavelength=%f4 m \nWaveNumber=%f2\n',...
    c,f0,w,T,lambda,k); 



