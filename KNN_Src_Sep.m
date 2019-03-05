%% Load the Data
load('eeg.mat')
first_channel=x_train(:,1,:);
second_channel=x_train(:,2,:);
third_channel=x_train(:,3,:);

S=size(first_channel);
f=reshape(first_channel,[S(1)*S(2),S(3)]);
s=reshape(second_channel,[S(1)*S(2),S(3)]);
t=reshape(third_channel,[S(1)*S(2),S(3)]);
%% Perform STFT for every channel and store the spectrogram of each 

%DFT length
N=64;

%Window size
nw=64;

%Hop size
h=48;

%Number of Frames will be same for every channel
p=1+floor((length(f)-nw)/h);

%Getting the Blackman window
w=blackman(nw);


%total samples on half range
ts=ceil((1+N)/2);

%Peforming STFT for every trial and appending data of these trials 
Xmk_first=stft(f(:,1),w,p,h,ts,nw,N);


%first channel
f1=abs(Xmk_first(3:7,:));
stft_size=size(f1);
f1=reshape(f1,[stft_size(1)*stft_size(2),1]);
final_f1=f1;

%second channel
Xmk_second=stft(s(:,1),w,p,h,ts,nw,N);
s1=abs(Xmk_second(3:7,:));
s1=reshape(s1,[stft_size(1)*stft_size(2),1]);
final_s1=s1;

%third channel
Xmk_third=stft(t(:,1),w,p,h,ts,nw,N);
t1=abs(Xmk_third(3:7,:));
t1=reshape(t1,[stft_size(1)*stft_size(2),1]);
final_t1=t1;

for i=2:size(f,2)
    
    Xmk_first=stft(f(:,i),w,p,h,ts,nw,N);
    
    f1=abs(Xmk_first(3:7,:));
    f1=reshape(f1,[stft_size(1)*stft_size(2),1]);
    final_f1=[final_f1 f1];
    
    Xmk_second=stft(s(:,i),w,p,h,ts,nw,N);
    s1=abs(Xmk_second(3:7,:));
    s1=reshape(s1,[stft_size(1)*stft_size(2),1]);
    final_s1=[final_s1 s1]
    
    
    Xmk_third=stft(t(:,i),w,p,h,ts,nw,N);
    t1=abs(Xmk_third(3:7,:));
    t1=reshape(t1,[stft_size(1)*stft_size(2),1]);
    final_t1=[final_t1 t1];

    
       
end 



%% Form data Matrix and Randomly initialize Projection Matrix
Z=[final_f1;final_s1;final_t1];
%projection matrix A
c=-1;
d=1;
%L dimensions
L=[2,4,6,8,10,15,20];
%K knn
K=[3,5,7,9,11,15,21];

%% Prepare the test dataset by taking stft of it 
test_first_channel=x_te(:,1,:);
test_second_channel=x_te(:,2,:);
test_third_channel=x_te(:,3,:);


S_te=size(test_first_channel);
fte=reshape(test_first_channel,[S_te(1)*S_te(2),S_te(3)]);
ste=reshape(test_second_channel,[S_te(1)*S_te(2),S_te(3)]);
tte=reshape(test_third_channel,[S_te(1)*S_te(2),S_te(3)]);

% Perform STFT for every channel and get the final test data matrix

%DFT length
N=64;

%Window size
nw=64;

%Hop size
h=48;

%Number of Frames will be same for every channel
p=1+floor((length(fte)-nw)/h);

%Getting the Blackman window
w=blackman(nw);


%total samples on half range
ts=ceil((1+N)/2);

%Peforming STFT for every trial and appending data of these trials 
Xmk_first_test=stft(fte(:,1),w,p,h,ts,nw,N);


%first channel
f1_test=abs(Xmk_first_test(3:7,:));
stft_size=size(f1_test);
f1_test=reshape(f1_test,[stft_size(1)*stft_size(2),1]);
final_f1_test=f1_test;

%second channel
Xmk_second_test=stft(ste(:,1),w,p,h,ts,nw,N);
s1_test=abs(Xmk_second_test(3:7,:));
s1_test=reshape(s1_test,[stft_size(1)*stft_size(2),1]);
final_s1_test=s1_test;

%third channel
Xmk_third_test=stft(tte(:,1),w,p,h,ts,nw,N);
t1_test=abs(Xmk_third_test(3:7,:));
t1_test=reshape(t1_test,[stft_size(1)*stft_size(2),1]);
final_t1_test=t1_test;

for i=2:size(fte,2)
    
    Xmk_first_test=stft(fte(:,i),w,p,h,ts,nw,N);
    
    f1_test=abs(Xmk_first_test(3:7,:));
    f1_test=reshape(f1_test,[stft_size(1)*stft_size(2),1]);
    final_f1_test=[final_f1_test f1_test];
    
    Xmk_second_test=stft(ste(:,i),w,p,h,ts,nw,N);
    s1_test=abs(Xmk_second_test(3:7,:));
    s1_test=reshape(s1_test,[stft_size(1)*stft_size(2),1]);
    final_s1_test=[final_s1_test s1_test];
    
    
    Xmk_third_test=stft(tte(:,i),w,p,h,ts,nw,N);
    t1_test=abs(Xmk_third_test(3:7,:));
    t1_test=reshape(t1_test,[stft_size(1)*stft_size(2),1]);
    final_t1_test=[final_t1_test t1_test];

    
       
end 

Ztest=[final_f1_test;final_s1_test;final_t1_test];

%%
accuracy=zeros(size(L,1),size(K,1));
for lval=1:length(L)
    for kval=1:length(K)
        A=c+(d-c).*rand(L(lval),size(Z,1));
        A=normalize(A);
        Y=sign(A*Z);
        Y(Y==(-1))=0;
        
        Atest=c+(d-c).*rand(L(lval),size(Ztest,1));
        Atest=normalize(Atest);
        Ytest=sign(Atest*Ztest);
        Ytest(Ytest==(-1))=0;
        
        y_pred=perform_knn(Ytest,Y,y_train,K(kval));
        accuracy(lval,kval)=get_accuracy(y_pred,y_te);
    end
    
end

disp(accuracy)

%% Perform Knn 
function y_pred=perform_knn(Ytest,Y,y_train,k)

Y=transpose(Y);
Ytest=transpose(Ytest);
dis=pdist2(Ytest,Y,'hamming');
[dis,idx]=sort(dis,2,'ascend');
dis=dis(:,1:k);
idx=idx(:,1:k);
y_pred=mode(y_train(idx),2);


end

%% Accuracy of the data 
function accuracy=get_accuracy(y_pred,y_te)
    accuracy=sum(y_pred==y_te)/length(y_te);
end

%% Normalize A so that every row is unit vector of it
function H=normalize(A)
for i=1:size(A,1)
    H(i,:)=A(i,:)./norm(A(i,:));
end
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











