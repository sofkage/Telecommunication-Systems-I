%Sofia Kafritsa Georganta 2016030136

clear all; 
close all;
clc;

%--------------Erwthma A.1-------------------
ka = 1/2*3

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

bits=(sign(randn(N,1))+1)/2;
Xn = bits_to_2PAM(bits);
X_n=1/Ts*upsample(Xn,over);
t2 = 0:Ts:N*T-Ts;
Xt =conv(X_n,phi)*Ts;
t_conv = min(t1)+min(t2):Ts:max(t1)+max(t2);

thita = 2*pi*rand(1);
fo = 2000;
Y_t = Xt.*( cos(2*pi*fo*t_conv + thita));

T_total=length(Y_t)*Ts; %
YF=fftshift(fft(Y_t,Nf))*Ts;  %fourier
PyF=(abs(YF).^2)/T_total; %mias ylopoihshs

PyF_pinakas = zeros(1000,Nf);

for i=1:1000
    
bits=(sign(randn(N,1))+1)/2;
Xn = bits_to_2PAM(bits);
X_n=1/Ts*upsample(Xn,over);
t2 = 0:Ts:N*T-Ts;
Xt =conv(X_n,phi)*Ts;
t_conv = min(t1)+min(t2):Ts:max(t1)+max(t2);

thita = 2*pi*rand(1);

Y_t = Xt.* cos(2*pi*fo*t_conv + thita);

YF=fftshift(fft(Y_t,Nf))*Ts;  %fourier
PyF=(abs(YF).^2)/T_total; %mias ylopoihshs

PyF_pinakas(i,:)=(abs(YF).^2)/T_total; %mias ylopoihshs

end

Sy_est = sum(PyF_pinakas)./1000;  %estimated

%var(Xn)=1
Rx = 1/T*conv(phi,phi)*Ts;
t_b = -2*A*T:Ts:2*A*T;
Ry = (1/2)*cos(2*pi*fo*t_b).*Rx;

Sy_theor = abs(fftshift(fft(Ry,Nf)*Ts));

%plot 
figure();
plot(f,PyF)
title('Plot of one periodogram with 2-PAM configuration');
xlabel('F');
ylabel('Py(F)')
grid on;

   %plot 
figure();
semilogy(f,PyF)
title('Semilogy of one periodogram with 2-PAM configuration');
xlabel('F');
ylabel('Py(F)')
grid on;


figure();
semilogy(f,Sy_est,'r')
title('Estimated PSD with 2-PAM configuration');
xlabel('F');
ylabel('Sy(F)');
legend('Sy,estimated')
grid on;

figure();
semilogy(f,Sy_est,'r',f,Sy_theor,'g')
title('Theoretical and Estimated PSD with 2-PAM configuration');
xlabel('F');
ylabel('Sy(F)');
legend('Sy,estimated', 'Sy, theoretical')
grid on;