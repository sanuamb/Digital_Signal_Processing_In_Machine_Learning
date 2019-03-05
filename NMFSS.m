%Single Channel Source Separation using NMF

%% Training Speech Signal

%Read the signal
[sp,fsp]=audioread('trs.wav');

%%DFT and window length
N=1024;
nw=1024;

%Hop size
h=N/2;

%%Number of frames
p=1+floor((length(sp)-nw)/h);

%Getting the hann window
w=hann(nw);

%total samples on half range
ts=ceil((1+N)/2);

Xmks=stft(sp,w,p,h,ts,nw,N);

%Calculate the time and Frequency vector
%the time vector for the spectrum is difference between the window length
%and hop size
%t=(nw/2:h:nw/2+(p-1)*h)/fs;

%frequency vector
%f=(0:ts-1)*fs/N;

%plotting the power spectogram
%figure;
%imagesc(t,transpose(f),(abs(Xmk)));
%title('Spectogram (db) of Speech Wave')
%xlabel('Time (sec)')
%ylabel('Frequency (Hz)');

%Taking the magnitudes S
S=abs(Xmks);
%% Get the weights of Speech Signal 

%Initialize W and H with random values
W=rand(size(Xmks,1),30);

H=rand(30,size(Xmks,2));

%Training NMF function
Ws=train_NMF(S,W,H);

%% Train The Noise Signal

%Read the signal
[nn,fsn]=audioread('trn.wav');

%%DFT and window length
N=1024;
nw=1024;

%Hop size
h=N/2;

%%Number of frames
p=1+floor((length(nn)-nw)/h);

%Getting the hann window
w=hann(nw);

%total samples on half range
ts=ceil((1+N)/2);

Xmk=stft(nn,w,p,h,ts,nw,N);

%Calculate the time and Frequency vector
%the time vector for the spectrum is difference between the window length
%and hop size
%t=(nw/2:h:nw/2+(p-1)*h)/fs;

%frequency vector
%f=(0:ts-1)*fs/N;

%plotting the power spectogram
%figure;
%imagesc(t,transpose(f),(abs(Xmk)));
%title('Spectogram (db) of Speech Wave')
%xlabel('Time (sec)')
%ylabel('Frequency (Hz)');

%Taking the magnitudes S
S=abs(Xmk);
%% Get the Weights of the Noise Signal

%Initialize W and H with random values
W=rand(size(Xmk,1),30);

H=rand(30,size(Xmk,2));

%Training NMF function
Wn=train_NMF(S,W,H);

%% Speech+Noise Signal

%read the audio
[sn,fssn]=audioread('x_nmf.wav');

%%DFT and window length
N=1024;
nw=1024;

%Hop size
h=N/2;

%%Number of frames
p=1+floor((length(sn)-nw)/h);

%Getting the hann window
w=hann(nw);

%total samples on half range
ts=ceil((1+N)/2);

%getting the spectrogram
Xmksn=stft(sn,w,p,h,ts,nw,N);

%Magnitude of the spectrogram
Y=abs(Xmksn);

%% Combine the Weights Obtained above of Speech and Noise

W=[Ws Wn];

H=rand(60,size(Xmksn,2));

Hsn=trainSN_NMF(Y,W,H);

num=Ws*Hsn(1:30,:);
denom=(Ws*Hsn(1:30,:))+(Wn*Hsn(31:60,:));
Mk=num./denom;

%Recover Speech Source
Shat=Mk.*Xmksn;

%% Peforming ISTFT

%Plotting Spectrogram
%Calculate the time and Frequency vector
%the time vector for the spectrum is difference between the window length
%and hop size
t=(nw/2:h:nw/2+(p-1)*h)/fssn;

%frequency vector
f=(0:ts-1)*fssn/N;

%plotting the power spectogram
figure;
imagesc(t,transpose(f),(abs(Xmksn)));
title('Spectogram (db) of Speech Wave')
xlabel('Time (sec)')
ylabel('Frequency (Hz)');


xt=ISTFT(Shat,w,h,nw,N,p,length(sn));

%Displaying time domain representation
%Time domain representation
figure;
subplot(2,1,1);
Nlen=length(sn);
ls=Nlen/fssn;
t1=linspace(0,ls,Nlen);

plot(t1,sn);
title('Actual Time Domain Signal');
xlabel('Time')
ylabel('Speech');



%plotting reconstructed signal
%Time domain representation
subplot(2,1,2);
Rlen=length(xt);
lr=Rlen/fssn;
t2=linspace(0,lr,Rlen);

plot(t2,xt);
title('Reconstructed Time Domain Signal with suppressed noise');
xlabel('Time')
ylabel('Speech');

%Reconstruct the signal
audiowrite('reconNMF.wav',xt,fssn)

%% NMF for Speech+Noise Signal to update only Basis

function H=trainSN_NMF(S,W,H)
cnt=0;
Xhat=W*H;
Xhat=Xhat+eps;
div=S./Xhat;
ldiv=log(div);
l1=S.*ldiv;
lfinal=l1-S+Xhat;

while 1
    disp(cnt);
    
    old_sum=sum(sum(lfinal,2));
    
    %update
    Hnew=update_H(S,W,H);
    H=Hnew;
    
    %Calculate the error
    Xhat=W*H;
    Xhat=Xhat+eps;
    div=S./Xhat;
    ldiv=log(div);
    l1=S.*ldiv;
    lfinal=l1-S+Xhat;

    
    disp(sum(sum(lfinal,2)));
    %Check for convergence
    if old_sum-(sum(sum(lfinal,2)))<1
        break
    end
       
    cnt=cnt+1;
    
end

end



%% NMF 

function W=train_NMF(S,W,H)
%threshold=10000;
%obj fn: Xij.log(Xij)/Xhatij-Xij+Xhatij
%Xhat=W*transpose(H);
cnt=0;
Xhat=W*H;
Xhat=Xhat+eps;
div=S./Xhat;
ldiv=log(div);
l1=S.*ldiv;
lfinal=l1-S+Xhat;

while 1
    disp(cnt);
    
    old_sum=sum(sum(lfinal,2));
    
    %update
    Wnew=update_W(S,W,H);
    Hnew=update_H(S,Wnew,H);
    W=Wnew;
    H=Hnew;
    
    %Calculate the error
    Xhat=W*H;
    Xhat=Xhat+eps;
    div=S./Xhat;
    ldiv=log(div);
    l1=S.*ldiv;
    lfinal=l1-S+Xhat;

    
    disp(sum(sum(lfinal,2)));
    %Check for convergence
    if old_sum-(sum(sum(lfinal,2)))<1
        break
    end
       
    cnt=cnt+1;
    
end

end
%% Update Weights
function Wfinal=update_W(S,W,H)
    xd=W*H;
    o=ones(size(S,1),size(S,2));
    xd=xd+eps;
    j1=S./xd;
    num=j1*transpose(H);
    denom=o*transpose(H);
    denom=denom+eps;
    d1=num./denom;
    Wfinal=W.*d1;
 
end
%% Update Basis
function Hfinal=update_H(S,W,H)
    xd=W*H;
    o=ones(size(S,1),size(S,2));
    xd=xd+eps;
    j1=S./xd;
    num=transpose(W)*j1;
    denom=transpose(W)*o;
    denom=denom+eps;
    d1=num./denom;
    Hfinal=H.*d1;
end




%% STFT Procedure

%stft function
function xwd=stft(x,w,p,h,ts,nw,N)
%initialize a ptr to increase the hop size
ptr=0;

%Zero vector for padding
z=zeros(N,1);


for m=1:p
    %extract frame of Input Data
    xt=x(ptr+1:ptr+nw);
    %apply window to the input frame
    xw=xt.*w;
    
    %Zero padding since inbuilt fft function does it.
    %size(xw)
    xw(numel(z))=0;
    
    %Apply dft to the windowed frame
    df=mydft(xw);
    
    %taking only 512 samples 
    xwd(:,m)=df(1:ts);
    
    %Advance the counter by the hopsize
    ptr=ptr+h;
    
end


end



%dft function

function xn=mydft(xw)
    
N1=length(xw);
xn=zeros(1,N1);
nn=zeros(1,N1);

for i=1:N1
    for j=1:N1
        nn(j)=xw(j)*exp(-1i*2*pi*(i-1)*(j-1)/N1);
    end
    xn(i)=sum(nn);
end  
end



%% Performing the ISTFT

function xt=ISTFT(Xpk,w,l,nw,N,p,len)

xt=zeros(1,len);
%intialize the ptr
cntr=0;

for m=1:p
    
    %Get the input sequence frame by frame
    xx=Xpk(:,m);
    
    %Symmetric conjugate 
    xx= [xx; conj(xx(end-1:-1:2))];
    
    %Take the IDFT of the input
    xidft=real(myIDFT(xx));
    xidft=xidft(1:nw);
    
    %Getting the actual signal (Add)
    xt(cntr+1:cntr+nw)=xt(cntr+1:cntr+nw)+xidft;
    cntr=cntr+l;
    
end

end



%IDFT
function xt=myIDFT(Xpk)

N1=length(Xpk);

xt=zeros(1,N1);
nn=zeros(1,N1);

for i=1:N1
    for j=1:N1
        nn(j)=Xpk(j)*exp(-1i*2*pi*(-(i-1))*(j-1)/N1);
    end
    xt(i)=sum(nn)/N1;
end

end


