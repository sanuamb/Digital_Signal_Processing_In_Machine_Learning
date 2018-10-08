
%% 
[vs,fss]=audioread('speech.wav');
[vn,fsn]=audioread('noise.wav');

%% Question 3
%--------------Question 3----------------------------

%Part 1 and 2
Xks=myDFT(vs);
Xkn=myDFT(vn);

%part 3
%Displaying the magnitude and phase response of the Speech and Noise wave
magphase(Xks,Xkn,fss,fsn);
%--------------------------------------------------------------------

%% Question 4
%-------------Question 4------------------------------------------

%part 1
%Adding speech and Noise DFT (Frequency Domain)
Xkns=Xks+Xkn;

%Displaying the magnitude response of speech, noise, speech+noise
magresponse(Xks,Xkn,Xkns,fss,fsn);

%part 2

%Displaying time domain Representation of Speech, Noise and Speech+Noise
vsn=vs+vn;
TDRepresentation(vs,vn,vsn,fss,fsn)

%part 3
%Computing the DFT of Speech + Noise for time domain representation
Xknst=myDFT(vsn);

%Displaying the DFTs for the speech+noise in Time domain, speech+Noise in
%frequency Domain and DFT of Difference between time domain and Frequency
%Domain speech+Noise
DFTSpeechNoise(Xkns,Xknst,fss)

%----------------------------------------------------------------------
%% Question 5
%-------------Question 5---------------------------------------

%part 1 and 2
%IDFT

xts=myIDFT(Xks);
xtn=myIDFT(Xkn);

%part 3
%Display the original time domain signals and the IDFTs of the respective
DisplayIDFT(xts,xtn,vs,vn,fss,fsn)

%part 4
%Play the signals of Speech and Nouse 

%Original speech
soundsc(vs,fss)

%Reconstructed speech
sound(real(xts),fss);

%Original Noise
soundsc(vn,fsn)

%Reconstructed Noise
soundsc(real(xtn),fsn);
%%
function Xk=myDFT(xn)
N=length(xn);
Xk=zeros(1,N);
enk=zeros(1,N);

for j=1:N
    for i=1:N
        enk(i)=xn(i)* exp(-1i*2*pi*(i-1)*(j-1)/N);
    end
   Xk(j)=sum(enk);
end

end

function TDRepresentation(vs,vn,vsn,fss,fsn)

figure;
%speech length in secs and time vector
Ns=length(vs);
lss=Ns/fss;
ts=linspace(0,lss,Ns);

%Noise length in secs and time vector
Nn=length(vn);
lns=Nn/fsn;
tn=linspace(0,lns,Nn);

% Speech + Noise length and time vector
Nsn= length(vsn);
%both the audio signals have same sampling frequency
lsns=Nsn/fss;
tsn=linspace(0,lsns,Nsn);

subplot(3,1,1)
plot(ts,vs)
title('Time Domain Representation For Speech')
xlabel('Time')
ylabel('Speech')

subplot(3,1,2)
plot(tn,vn)
title('Time Domain Representation For Noise')
xlabel('Time')
ylabel('Noise')

subplot(3,1,3)
plot(tsn,vsn)
title('Time Domain Representation For Speech+Noise')
xlabel('Time')
ylabel('Speech+Noise')

end

function DFTSpeechNoise(Xkns,Xknst,fss)

figure;
%plotting the magnitude response of the Noise+Speech from the frequency
%Domain
subplot(3,1,1);
mag_snf=10*log10(abs(Xkns));
fsnf=(0:length(Xkns)-1)*fss/length(Xkns);
plot(fsnf,mag_snf);
title('DFT of Speech+Noise Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude (db)');

%plotting the magnitude response of the noise+speech from the time domain
subplot(3,1,2)
mag_snt=10*log10(abs(Xknst));
fsnt=(0:length(Xknst)-1)*fss/length(Xknst);
plot(fsnt,mag_snt);
title('DFT of Speech+Noise Time Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude (db)');

%plotting the magnitude response for the difference between the time domain
%and frequency domain DFTs
subplot(3,1,3);
a1=abs(Xknst);
a2=abs(Xkns);
Diff=(a1-a2);
mag_snd=((Diff));
fd=(0:length(Diff)-1)*fss/length(Diff);
plot(fd,mag_snd);
title('DFT of Difference between time Domain and Freq Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude (db)');


end


function magphase(Xks,Xkn,fss,fsn)

figure;
%Magnitude Response of Speech wave
subplot(2,2,1)
mag_speech=10*log10(abs(Xks));
fs=(0:length(Xks)-1)*fss/length(Xks);
plot(fs,mag_speech);
title('Magnitude response of Speech Wave')
xlabel('Frequency (Hz)')
ylabel('Magnitude (DB)');

%Magnitude Response of Noise wave
subplot(2,2,2)
mag_Noise=10*log10(abs(Xkn));
fn=(0:length(Xkn)-1)*fsn/length(Xkn);
plot(fn,mag_Noise);
title('Magnitude response of Noise Wave')
xlabel('Frequency (Hz)')
ylabel('Magnitude (DB)');

%Phase response of Speech wave
subplot(2,2,3)
phase_speech=angle(Xks);
plot(fs,phase_speech);
title('Phase response of Speech Wave')
xlabel('Frequency (Hz)')
ylabel('Phase');


%Phase response of Noise wave
subplot(2,2,4)
phase_noise=angle(Xkn);
plot(fn,phase_noise);
title('Phase response of Noise Wave')
xlabel('Frequency (Hz)')
ylabel('Phase');
end

function xt=myIDFT(Xkk)

N=length(Xkk);

xt=zeros(1,N);
nn=zeros(1,N);

for i=1:N
    for j=1:N
        nn(j)=Xkk(j)*exp(-1i*2*pi*(-(i-1))*(j-1)/N);
    end
    xt(i)=sum(nn)/N;
end

end

function DisplayIDFT(xts,xtn,vs,vn,fss,fsn)
%speech length in secs and time vector
Ns=length(vs);
lss=Ns/fss;
ts=linspace(0,lss,Ns);

%Noise length in secs and time vector
Nn=length(vn);
lns=Nn/fsn;
tn=linspace(0,lns,Nn);

%reconstructed length in secs and time vector of Speech
Nrs=length(xts);
lrs=Nrs/fss;
trs=linspace(0,lrs,Nrs);

%reconstructed length in secs and time vector of Noise
Nrn=length(xtn);
lrn=Nrn/fsn;
trn=linspace(0,lrn,Nrn);

figure;
subplot(2,1,1)
plot(ts,vs)
title('Original Time Domain Representation For Speech')
xlabel('Time')
ylabel('Speech')

subplot(2,1,2)
plot(trs,real(xts))
title('Reconstructed Time Domain Representation For Speech')
xlabel('Time')
ylabel('Speech')



figure;
subplot(2,1,1)
plot(tn,vn)
title('Original Time Domain Representation For Noise')
xlabel('Time')
ylabel('Noise')

subplot(2,1,2)
plot(trn,real(xtn))
title('Reconstructed Time Domain Representation For Noise')
xlabel('Time')
ylabel('Noise')


end

%Question 4_1
function magresponse(Xks,Xkn,Xkns,fss,fsn)
figure;
%Magnitude Response of Speech wave
subplot(3,1,1)
mag_speech=10*log10(abs(Xks));
fs=(0:length(Xks)-1)*fss/length(Xks);
plot(fs,mag_speech);
title('Magnitude response of Speech Wave (Frequency Domain)')
xlabel('Frequency (Hz)')
ylabel('Magnitude (DB)');

%Magnitude Response of Noise wave
subplot(3,1,2)
mag_Noise=10*log10(abs(Xkn));
fn=(0:length(Xkn)-1)*fsn/length(Xkn);
plot(fn,mag_Noise);
title('Magnitude response of Noise Wave (Frequency Domain)')
xlabel('Frequency (Hz)')
ylabel('Magnitude (DB)');

%Magnitude Response of Speech+Noise wave
subplot(3,1,3);
mag_sn=10*log10(abs(Xkns));
fns=(0:length(Xkns)-1)*fsn/length(Xkns);
plot(fns,mag_sn);
title('Magnitude response of Speech+Noise Wave (Frequency Domain)')
xlabel('Frequency (Hz)')
ylabel('Magnitude (DB)');

end