%Sofia Kafritsa Georganta 2016030136

clear all; 
close all;
clc;

%--------------Erwthma A.1-------------------

%orismos parametrwn
T = 10 ^ (-2);
over = 10; 
Ts = T/over;
A = 4;
a=[0 0.5 1];

%SRRC pulses me vash th synarthsh srrc_pulse

[phi1, t1] = srrc_pulse (T, Ts, A, a(1));   %gia a=0
[phi2, t2] = srrc_pulse (T, Ts, A, a(2));   %gia a=0.5
[phi3, t3] = srrc_pulse (T, Ts, A, a(3));   %gia a=1

%koino plot olwn twn palmwn
figure(1);
plot(t1,phi1,'b')
hold on
plot(t2,phi2,'g')
hold on
plot(t3,phi3,'r')
grid on
legend('a = 0','a = 0.5','a = 1');
hold off

title('SRRC Pulses ')
xlabel('Time in sec')
ylabel('SRRC Amplitude')


%--------------Erwthma A.2-------------------

Fs=1/Ts;
Nf=2048;

%orizw ton aksona twn syxnothtwn 
F = ((-Fs/2):(Fs/Nf):(Fs/2 - Fs/Nf));

%orizw tous metasxhmatismous fourier
fourier_1 = fftshift(fft(phi1,length(F))*Ts);
fourier_2 = fftshift(fft(phi2,length(F))*Ts);
fourier_3 = fftshift(fft(phi3,length(F))*Ts);

%orizw thn  fasmatikh pyknothta energeias  |PHI()|^2

ESD_1=abs(fourier_1).^2;
ESD_2=abs(fourier_2).^2;
ESD_3=abs(fourier_3).^2;

%koino plot fasmatikwn pyknothtwn energeias 
figure(2);
plot(F,ESD_1,'b')
hold on
plot(F,ESD_2,'g')
hold on
plot(F,ESD_3,'r')
legend('a = 0','a = 0.5','a = 1');
grid on;
hold off

title('Energy Spectral Density Plot')
xlabel('Frequency in Hz')
ylabel('EDS')

%koino semilogy twn esd kai diaoretikos tropos dhlwshs xwris hold on
figure(3);
semilogy(F,ESD_1,'b',F,ESD_2,'g',F,ESD_3,'r')
legend('a = 0','a = 0.5','a = 1');
grid on;

title('Energy Spectral Density Semilogy')
xlabel('Frequency in Hz')
ylabel('EDS')



%--------------Erwthma A.3-------------------

%ypologizw mesa se or loop gia to thewrhtiko bw gia tous 3 palmous

for a=0:0.5:1
  bw_theory = (1+a)/(2*T)
end

%orizw tis orizonties grammes 
c1 = T/(10^3);
c2 = T/(10^5);

%dhmiourgw tis grammes
line = ones (1,length(F));
line1=line*c1;
line2=line*c2;

figure(4);
semilogy(F,ESD_1,'b')
hold on
semilogy(F,ESD_2,'g')
hold on
semilogy(F,ESD_3,'r')
hold on
semilogy(F,line1)
hold on
semilogy(F,line2)
legend('a = 0','a = 0.5','a = 1', 'c1=T/(10^3)', 'c2=T/(10^5)');
grid on;
hold off

title('EDS with c')
xlabel('Frequency')
ylabel('ESD')


%--------------Erwthma B-------------------

A=5;
a=[0 0.5 1];
%thetoume A=5 kai ftiaxnoyme 3 kainourious palmous phi
[phi1_b, t1] = srrc_pulse (T, Ts, A, a(1));     %gia a=0
[phi2_b, t2] = srrc_pulse (T, Ts, A, a(2));   %gia a=0.5
[phi3_b, t3] = srrc_pulse (T, Ts, A, a(3));     %gia a=1



for k=[0 1 2 4]
    
    %dhniourgoume tous kathisterhmenous palmous, me mhdenika prin kai meta
    %to phi

    del_filter1(k+1,:)=[zeros(1,k*over) phi1_b zeros(1,(2*A-k)*over)];
    del_filter2(k+1,:)=[zeros(1,k*over) phi2_b zeros(1,(2*A-k)*over)];
    del_filter3(k+1,:)=[zeros(1,k*over) phi3_b zeros(1,(2*A-k)*over)];

    %orizoume koino aksona xronou
    t_total=[-A*T:Ts:A*(1+2)*T];

    %ypologismos ginomenwn
    product_1 = del_filter1(1,:).*del_filter1(k+1,:); 
    product_2 = del_filter2(1,:).*del_filter2(k+1,:); 
    product_3 = del_filter3(1,:).*del_filter3(k+1,:); 
    
    %ypoloogismos oloklhrwmatwn
    integral_1(k+1)= sum(product_1*Ts) 
    integral_2(k+1)= sum(product_2*Ts) 
    integral_3(k+1)= sum(product_3*Ts) 


end

%sxediasmos phi
figure(5)
plot(t_total,del_filter1,'b')
grid on
title('phi(t) and phi(t-kT) for a=0 ');
xlabel('Time in sec');
ylabel('Pulses');

figure(6)
plot(t_total,del_filter2,'g')
grid on
title('phi(t) and phi(t-kT) for a=0.5 ');
xlabel('Time in sec');
ylabel('Pulses');

figure(7)
plot(t_total,del_filter3,'r')
grid on
title('phi(t) and phi(t-kT) for a=1 ');
xlabel('Time in sec');
ylabel('Pulses');

%sxediasmos ginomenwn

figure(8)
plot(t_total,product_1,'b')
grid on
title(['Product 1 for a=0 ']);
xlabel('Time in sec');
ylabel('Product');

figure(9)
plot(t_total,product_2,'g')
grid on
title(['Product 2 for a=0.5 ']);
xlabel('Time in sec');
ylabel('Product');


figure(10)
plot(t_total,product_3,'r')
grid on
title(['Product 3 for a=1 ']);
xlabel('Time in sec');
ylabel('Product');





%--------------Erwthma C.1-------------------

N=100;
T=0.1;
Ts=T/over;
a=0.5;
b=(sign(randn(N,1))+1)/2;

%--------------Erwthma C.2-------------------
%(a) 
X = bits_to_2PAM(b);

%(b)
X_delta=1/Ts*upsample(X,over);
t_delta = 0:Ts:N*T-Ts;

%Sxediash
figure
plot(t_delta,X_delta)
title('X delta Signal');
xlabel('Time in sec');
ylabel('X_d');
grid on

%(c)

%Dhmiourgoume ton palmo SRRC me ta nea dedomena
[phi_apokommenos, t_phi_apok] = srrc_pulse(T, Ts, A, a);

%ypologizoume th syneliksi tou Xdelta me ton apokommeno phi
X_t = Ts*conv(X_delta,phi_apokommenos); 

%orizoume ton aksona tou xronou
t_conv1 = (min(t_phi_apok)+min(t_delta):Ts:max(t_phi_apok)+max(t_delta));

figure
plot(t_conv1,X_t)
title('Convolution X(t)=Xdelta(t)*phiApokommenos(t)');
xlabel('Time in sec');
ylabel('X_conv');
grid on


%(d)
Z_t = Ts*conv(X_t,phi_apokommenos);
%t_conv2 = [-2*A*T:Ts:2*A*T+length([X])-Ts];
t_conv2 = (min(t_conv1)+min(t_phi_apok):Ts:max(t_conv1)+max(t_phi_apok));

figure
plot(t_conv2,Z_t)
hold on
stem((0:N-1)*T,X,'r');
title('Convolution Z(t)=X(t)*phiApokommenos(-t)');
xlabel('Time in sec');
ylabel('Z_conv');
legend('Z(t)', 'Xk')
grid on






