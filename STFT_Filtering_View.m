
%STFT using filtering View

%Read the input signal
[v,fs]=audioread('speech2.wav');

%Input the parameters from the user

%Length of DFT
prompt1='Enter the DFT length(N):';
N=input(prompt1);

%Window Length
prompt2='Enter the window Length:';
Nw=input(prompt2);

%Window Length
nw=ceil(Nw*fs); %Window length in samples


%Getting the Hamming window vector
w=hamming(nw);



%length of the signal
len=length(sp);


%Number of samples (half of N)
rn=ceil((1+N)/2);


Xpk=stft(sp,w,nw,N);



%Calculate the time and Frequency vector
%the time vector for the spectrum is difference between the window length
%and hop size

t=linspace(0,length(Xpk)/fs,length(Xpk));

%frequency vector
f=(0:rn-1)*fs/N;

%plotting the power spectogram
figure;
imagesc(t,transpose(f),20*log10(abs(Xpk)));
title('Spectogram (db) of Speech Wave')
xlabel('Time (sec)')
ylabel('Frequency (Hz)');



%STFT function
function xwd=stft(x,w,nw,N)

%exp1=(-j*2*pi)/rn;
%e=exp(exp1);
n=1:nw;
xwd=zeros(N,length(x)+nw-1);
exp1=exp(1i*2*pi/N);

for k=0:N-1
    e1=exp1.^(n.*k);
    
    w1=w.*e1';
    
    xcon=conv(x,w1);
    
    lxcon=1:length(xcon);
    
    e2=exp1.^(lxcon.*(-k));
    
    m=xcon.*e2';
    
    xwd(k+1,:)=m;   
end
xwd=xwd(1:end/2,:);
end


