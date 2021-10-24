%Sofia Kafritsa Georganta 2016030136

clear all; 
close all;
clc;

%--------------Erwthma A.1-------------------

%orismos parametrwn
T = 2*(10 ^ (-3));
over = 20; 
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

figure();
semilogy(f,esd);
grid on;
title('ESD');
xlabel('F(Hz)');
ylabel('|PHI(F)|^2');

%--------------Erwthma A.5-------------------

N=100;

b=(sign(randn(N,1))+1)/2;
Xn = bits_to_2PAM(b);
X_n=1/Ts*upsample(Xn,over);
t2 = 0:Ts:N*T-Ts;
Xt =conv(X_n,phi)*Ts;
t_conv = min(t1)+min(t2):Ts:max(t1)+max(t2);


figure()
plot(t_conv,Xt)
title('X(t) from 2-PAM');
xlabel('Time in sec');
ylabel('X(t)');
grid on



K=1000;
T_total=length(Xt)*Ts; %
XF=fftshift(fft(Xt,Nf))*Ts;  %fourier
PxF=(abs(XF).^2)/T_total %mias ylopoihshs


%ftiaxnoume polles ylopoihseis gia thn Sx
for k=1:K
    
    b=(sign(randn(N,1))+1)/2;
    Xn = bits_to_2PAM(b);
    X_n=1/Ts*upsample(Xn,over);
    Xt = conv(X_n,phi)*Ts;
    XF=fftshift(fft(Xt,Nf))*Ts;
    T_total=length(Xt)*Ts; %

    PxF=(abs(XF).^2)/T_total;
    PxF_pinakas(k,:)=PxF;

end


Sx_est = sum(PxF_pinakas)./K;  %estimated
Sx_theor = (1/T).*(esd); %theoritiki timi


%plot 
figure();
plot(f,PxF)
title('Plot of one periodogram with 2-PAM configuration and T=2T');
xlabel('F');
ylabel('Px(F)')
grid on;


%semilogy
figure();
semilogy(f, PxF);
title('Semilogy of periodogram with 2-PAM configuration and T=2T');
xlabel('F');
ylabel('Px(F)')
grid on;
    
    
figure();
semilogy(f,Sx_est,'r',f,Sx_theor,'g')
title('Theoretical and Estimated PSD with 2-PAM configuration');
xlabel('F');
ylabel('Sx(F)');
legend('Sx,estimated', 'Sx, theoretical')
grid on;



