clear all;
clc;

%% Question 1

%Single Glottal Pulse

%given sampling rate
fs=10000;

%vector n
n=-0.05:1/fs:0.05;

%alpha
alpha=0.99;

%getting the result from the given equation (time reversed decaying exponential)

%pre-allocating memory
u=zeros(1,length(n));
res=getresult(u,n,alpha);

%Plotting the result as function of n
figure;
plot(n,res);
title('Time reversed Decaying exponential');
xlabel('Time(sec)')
ylabel('Amplitude')

%Performing convolution
gpulse=conv(res,res);

%Plotting the convolution results
figure;
n1=-0.1:1/fs:0.1;
plot(n1,gpulse);
title('Single glottal pulse');
xlabel('Time (sec)');
ylabel('Glottal Pusle');

%1.2 Glottal Pulse Train
%Sampling rate is already defined

%define vector t
t=-1:1/fs:1;

%frequency of impulse train
ft=10;

%Computing value for pitch period P
p=1/ft;

%setting the sample pulse train at t=-1 seconds to 1 and every pth sample
%thereafter to 1

%pre-allocating for faster execution
imp=zeros(1,length(t));

for j=1:p*fs:length(t)
    imp(j)=1;
end

%Plotting the signal using the stem function
figure;
stem(t,imp);
title('Impulse train');
xlabel('Time(sec)');
ylabel('Impulse');

%Convolve glottal pulse with the pusle train

res1=conv(gpulse,imp);

%Plotting the glottal Pulse Train
t1=-1.1:1/fs:1.1;
figure;
plot(t1,res1);
title('Convolved Glottal Pulse with impulse Train');
xlabel('Time(sec)');
ylabel('Pulse Amplitude');

%% Question 2 Simulating vocal tract model for a vowel

%given parameters as follows
formants=[530,1840,2480,4000];
bandwidths=[50,80,100,150];
fs=10000;
nfft=2048;

%creating impulse response
impres=genVowelImpulseResponse(formants,bandwidths,nfft,fs);

%plot the first 200 samples
len=length(impres);
Ns=len/fs;
ts=linspace(0,Ns,len);

figure;
plot(ts(1:200),impres(1:200));
title('Impulse response as a function of time');
xlabel('Time(sec)');
ylabel('Impulse');

%2.2 plot first half of the log-power magnitude spectrum of the impulse
%response
imp1=fft(impres,nfft);
logpow=10*log10(abs(imp1));
fp=(0:length(imp1)-1)*fs/length(imp1);

%plotting
figure;
plot(fp(1:length(fp)/2),logpow(1:length(imp1)/2));
title('log-power magnitude spectrum');
xlabel('Frequency (Hz)');
ylabel('log-power magnitude(db)');

%2.3 generate single glottal pulse using alternative approach
gp=genGlottalPulse(fs);
len1=length(gp);
Ng=len1/fs;
tgp=linspace(0,Ng,len1);

figure;
plot(tgp,gp);
title('Single glottal pulse with alternative approach');
xlabel('Time(sec)');
ylabel('Glottal Pulse');

%Generate pulse train

%pre-allocating
pt=zeros(1,1000);

for j=1:50:length(pt)
    pt(j)=1;
end

%Convolving the impulse train with the glottal pulse
res2=conv(pt,gp);

%plotting the result
tr=linspace(0,length(res2)/fs,length(res2));
figure;
plot(tr,res2);
title('Convolved Glottal pulse and impulse train');
xlabel('Time(sec)');
ylabel('Pulse Amplitude');

%Convolve previous result with vocal tract impulse response
res3=conv(res2,impres);
tgi=linspace(0,length(res3)/fs,length(res3));
figure;
plot(tgi,res3);
title('Convolved previous result with vocal tract impulse response');
xlabel('Time(sec)');
ylabel('Pulse Amplitude');

%play the result
soundsc(res3,fs);

%% Question 3 Gammatone filterbank

frange=[50 8000];
numFilters=32;
fs1=16000;
l=4;
N=2048;

erb_b=hz2erb(frange);
erb=[erb_b(1):diff(erb_b)/(numFilters-1):erb_b(2)];
f=erb2hz(erb);
b=1.019*24.7*(4.37*f/1000+1);

%given the equation for the gammatone filterbank generating it
vt=linspace(0,N/fs1,N);

%pre-allocating g
g=zeros(numFilters,N);

%generating the filterbanks
for i=1:numFilters
    g(i,:)=vt.^(l-1).* exp(-2 * pi * b(i)* vt).*cos(2 * pi * f(i)* vt);
end

%plotting the filterbanks 
%fg1=(0:N-1)*fs1/N;

%half of the frequency range for display 
fs2=0:fs1/N:fs1/2 - fs1/N;

figure;
for i=1:numFilters
    g1(i,:)=fft(g(i,:),N);
    mag=10*log10(abs(g1(i,1:N/2)));
    plot(fs2,mag);
    title('Gammatone Filterbank');
    xlabel('freq(Hz)');
    ylabel('Mag Gain (db)');
    hold on;
end


%Filtering
[s,fsp]=audioread('speech.wav');

%pre-allocating
ft=zeros(length(s),numFilters);

for i=1:numFilters
    ft(:,i)=fftfilt(g(i,:),s);
end

%plotting the speech signal
tsp=linspace(0,length(s)/fsp,length(s));

figure;

subplot(2,3,1);
plot(tsp,s);
title('Original Speech Signal');
xlabel('Time(sec)');
ylabel('Speech');

subplot(2,3,2);
plot(tsp,ft(:,1));
title('filtered by 1st filter');
xlabel('Time(sec)');
ylabel('Speech');

subplot(2,3,3);
plot(tsp,ft(:,8));
title('filtered by 8th filter');
xlabel('Time(sec)');
ylabel('Speech');

subplot(2,3,4);
plot(tsp,ft(:,16));
title('filtered by 16th filter');
xlabel('Time(sec)');
ylabel('Speech');

subplot(2,3,5);
plot(tsp,ft(:,24));
title('filtered by 24th filter');
xlabel('Time(sec)');
ylabel('Speech');

subplot(2,3,6);
plot(tsp,ft(:,32));
title('filtered by 32nd filter');
xlabel('Time(sec)');
ylabel('Speech');

%listen to the six signals plotted
soundsc(real(double(s)), fsp);
pause(5);
soundsc(real(double(ft(:,1))), fsp);
pause(5);
soundsc(real(double(ft(:,8))), fsp);
pause(5);
soundsc(real(double(ft(:,16))), fsp);
pause(5);
soundsc(real(double(ft(:,24))), fsp);
pause(5);
soundsc(real(double(ft(:,32))), fsp);
%% 
function x=getresult(u,n,alpha)

    x=zeros(1,length(n));

%from the given conditions 
    for i=1:length(n)
        if (n(i)<=0)
            u(i)=1;
        else
            u(i)=0;
        end
        x(i)=(alpha^-n(i))*u(i);
    end
end

    