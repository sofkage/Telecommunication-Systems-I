clear all; 
close all;
clc;

%--------------Erwthma A.1-------------------

%orismos parametrwn
T = 10 ^ (-3);
over = 10; 
Ts = T/over;
A = 4;
a = 0.5;

%SRRC pulse me vash th synarthsh srrc_pulse
[phi, t1] = srrc_pulse (T, Ts, A, a);   

Fs=1/Ts;
Nf=2048;

%aksonas syxnothtwn
f = (-Fs/2):(Fs/Nf):(Fs/2)-(Fs/Nf);

%fourier transform
fourier = fftshift(fft(phi,Nf))*Ts;

%esd
esd=abs(fourier).^2;
N=100;
K=1000;

%--------------Erwthma A.4-------------------

b=(sign(randn(N,1))+1)/2;
Xn = bits_to_4PAM(b);
X_n=1/Ts*upsample(Xn,over);
Xt =conv(X_n,phi)*Ts;

T_total=length(Xt)*Ts; %
XF=fftshift(fft(Xt,Nf))*Ts;  %fourier
PxF=(abs(XF).^2)/T_total %mias ylopoihshs

for k=1:K
    
    b=(sign(randn(N,1))+1)/2;
    Xn = bits_to_4PAM(b);
    X_n=1/Ts*upsample(Xn,over);
    Xt = conv(X_n,phi)*Ts;
    XF=fftshift(fft(Xt,Nf))*Ts;
    
    PxF=(abs(XF).^2)/T_total;
    PxF_pinakas(k,:)=PxF;

end


Sx_est = sum(PxF_pinakas)./K;  %estimated
Sx_theor = (var(Xn)/T).*(esd); %theoritiki timi

%plot 
figure();
plot(f,PxF)
title('Plot of one periodogram with 4-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;


%semilogy
figure();
semilogy(f, PxF);
title('Semilogy of periodogram with 4-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;
    
    
figure();
semilogy(f,Sx_est,'r',f,Sx_theor,'g')
title('Theoretical and Estimated PSD with 4-PAM configuration');
xlabel('F');
ylabel('Sx(F)');
legend('Sx,estimated', 'Sx, theoretical')
grid on;