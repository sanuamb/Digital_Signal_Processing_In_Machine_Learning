
%Implemented STFT Fourier View and ISTFT (Overlap and Add)

%Read the audio file
[sp,fs]=audioread('speech2.wav');

%Input the parameters from the user

%Length of DFT
prompt1='Enter the DFT length(N):';
N=input(prompt1);

%Window Length
prompt2='Enter the window Length:';
Nw=input(prompt2);

%Hop/Step Size
prompt3='Enter the Hop size:';
L=input(prompt3);

%Given Parameters will be in secs 

%Hop
l=floor(L*fs); %Step in samples

%Window Length
nw=ceil(Nw*fs); %Window length in samples


%Getting the Hamming window vector
w=hamming(nw);

%length of the signal
len=length(sp);

%Number of Frames p
p=1+floor((len-nw)/l);

%Number of samples (half of N)
rn=ceil((1+N)/2);

%STFT of the i/p signal
Xpk=stft(sp,w,l,nw,N,p,rn);

%ISTFT 
xt=ISTFT(Xpk,w,l,nw,N,p,len);


%Calculate the time and Frequency vector
%the time vector for the spectrum is difference between the window length
%and hop size
t=(nw/2:l:nw/2+(p-1)*l)/fs;



%frequency vector
f=(0:rn-1)*fs/N;

%plotting the power spectogram
figure;
imagesc(t,transpose(f),20*log10(abs(Xpk)));
title('Spectogram (db) of Speech Wave')
xlabel('Time (sec)')
ylabel('Frequency (Hz)');



%Displaying the actual vs reconstructed time signal 
%Actual
Nlen=length(sp);
ls=Nlen/fs;
ts=linspace(0,ls,Nlen);

%Reconstructed
Nlen1=length(xt);
lt=Nlen1/fs;
tt=linspace(0,lt,Nlen1);

%Plotting
figure;
subplot(2,1,1);
plot(ts,sp);
title('Actual Time Domain Signal');
xlabel('Time')
ylabel('Speech');

subplot(2,1,2);
plot(tt,xt);
title('Reconstructed Time Domain Signal');
xlabel('Time')
ylabel('Speech');



%STFT function
function xwd=stft(x,w,l,nw,N,p,rn)

%Intialize the pointer
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
    xwd(:,m)=df(1:rn);
    
    %Advance the counter by the hopsize
    ptr=ptr+l;
    
end

end


%DFT 
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


function xt=ISTFT(Xpk,w,l,nw,N,p,len)

xt=zeros(1,len);
%intialize the ptr
cntr=0;

for m=1:p
    
    %Get the input sequence frame by frame
    xx=Xpk(:,m);
    xx= [xx; conj(xx(end-1:-1:2))];
    
    %Take the IDFT of the input
    xidft=real(myIDFT(xx));
    xidft=xidft(1:nw);
    
    %Getting the actual signal (Add)
    xt(cntr+1:cntr+nw)=xt(cntr+1:cntr+nw)+xidft;
    cntr=cntr+l;
    
end

%Scaling the signal
W0=sum(w.^2);
xt=xt.*l/W0;
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



