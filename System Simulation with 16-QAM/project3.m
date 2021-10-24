clear all;
close all;

%Dedomena gia ta erwthmata
N = 200;
A = 1;
A_srrc = 4;
T = 0.01;
over = 10;
Ts = T/over;
a=0.5;

Fs=1/Ts;
Nf=2048;

%-----Erwthma 1-----  

%dhmiourgia duadikhs akolouthias me 4N isopithana bits
bit_seq = (sign(randn(4*N,1))+1)/2;

%-----Erwthma 2,3-----

XI=bits_to_4_PAM(bit_seq(1:2*N),A);  %ta prwta 2N bits
XQ=bits_to_4_PAM(bit_seq(2*N+1:4*N),A); %ta epomena 2N bits

%----Erwthma 4-----

[phi,t1]=srrc_pulse(T, Ts, A_srrc, a);

Xi_n = 1/Ts*upsample(XI,over);
Xi_t = conv(Xi_n,phi)*Ts;

Xq_n = 1/Ts*upsample(XQ,over);
Xq_t = conv(Xq_n,phi)*Ts;

t2 = 0:Ts:N*T-Ts;   
t_conv = min(t1)+min(t2):Ts:max(t1)+max(t2);

figure();
plot(t_conv,Xi_t);
grid on;
title('Output waveform Xi(t)');
xlabel('t(s)');
ylabel('Xi(t)');

figure();
plot(t_conv,Xq_t);
grid on;
title('Output waveform Xq(t)');
xlabel('t(s)');
ylabel('Xq(t)');


f = (-Fs/2):(Fs/Nf):(Fs/2)-(Fs/Nf);   

T_total_i=length(Xi_t)*Ts; 
XF_i=fftshift(fft(Xi_t,Nf))*Ts;  %fourier
PxF_i=(abs(XF_i).^2)/T_total_i %mias ylopoihshs

T_total_q=length(Xq_t)*Ts; 
XF_q=fftshift(fft(Xq_t,Nf))*Ts;  %fourier
PxF_q=(abs(XF_q).^2)/T_total_q %mias ylopoihshs

figure();
plot(f, PxF_i);
title('Plot of periodogram of Xi(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PxF_q);
title('Plot of periodogram of Xq(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

%------Erwthma 5---------

Fo = 300;
Xi_mod_t = Xi_t.*(2*cos(2*pi*Fo*t_conv))*Ts;
Xq_mod_t = Xq_t.*(-2*sin(2*pi*Fo*t_conv))*Ts;

figure();
plot(t_conv,Xi_mod_t);
grid on;
title('Xi(t) multiplied with 2cos(2piFot)');
xlabel('t(s)');
ylabel('Xi(t)');

figure();
plot(t_conv,Xq_mod_t);
grid on;
title('Xq(t) multiplied with -2sin(2piFot)');
xlabel('t(s)');
ylabel('Xq(t)');


T_total_i=length(Xi_mod_t)*Ts; %
XF_i=fftshift(fft(Xi_mod_t,Nf))*Ts;  %fourier
PxF_i=(abs(XF_i).^2)/T_total_i %mias ylopoihshs

T_total_q=length(Xq_mod_t)*Ts; %
XF_q=fftshift(fft(Xq_mod_t,Nf))*Ts;  %fourier
PxF_q=(abs(XF_q).^2)/T_total_q %mias ylopoihshs


figure();
plot(f, PxF_i);
title('Plot of periodogram of Xi(t)_m_o_d');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PxF_q);
title('Plot of periodogram of Xq(t)_m_o_d');
xlabel('F');
ylabel('Px(F)')
grid on;



%------Erwthma 6---------

X_t_mod = Xi_mod_t + Xq_mod_t;

figure();
plot(t_conv,X_t_mod);
grid on;
title('Input waveform X(t)');
xlabel('t(s)');
ylabel('X(t)');

T_total=length(X_t_mod)*Ts; %
XF=fftshift(fft(X_t_mod,Nf))*Ts;  %fourier
PxF=(abs(XF).^2)/T_total %mias ylopoihshs

figure();
plot(f, PxF);
title('Plot of periodogram of X(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

%-----Erwthma 8-----  

SNR = 30;
var_w = 10*(A^2)/(Ts*10^(SNR/10));
%var_n = Ts*var_w/2;

noise = sqrt(var_w)*randn(1,length(X_t_mod));

X_mod_noise = X_t_mod + noise;


%-----Erwthma 9-----  


Xi_demodulated = X_mod_noise.*cos(2*pi*Fo*t_conv)*Ts;
Xq_demodulated = X_mod_noise.*(-1*sin(2*pi*Fo*t_conv))*Ts;

figure();
plot(t_conv,Xi_demodulated);
grid on;
title('?i(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with cos(2piFot)');
xlabel('t(s)');
ylabel('Xi(t)');

figure();
plot(t_conv,Xq_demodulated);
grid on;
title('?q(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with -sin(2piFot)');
xlabel('t(s)');
ylabel('Xq(t)');

T_total_i=length(Xi_demodulated)*Ts; %
XF_i=fftshift(fft(Xi_demodulated,Nf))*Ts;  %fourier
PxF_i=(abs(XF_i).^2)/T_total_i %mias ylopoihshs

T_total_q=length(Xq_demodulated)*Ts; %
XF_q=fftshift(fft(Xq_demodulated,Nf))*Ts;  %fourier
PxF_q=(abs(XF_q).^2)/T_total_q %mias ylopoihshs


figure();
plot(f, PxF_i);
title('Plot of periodogram of Xi(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with cos(2piFot)');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PxF_q);
title('Plot of periodogram of ?q(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with -sin(2piFot)');
xlabel('F');
ylabel('Px(F)')
grid on;

%-----Erwthma 10-----  

Xi_demodulated = conv(Xi_demodulated,phi)*Ts;
Xq_demodulated = conv(Xq_demodulated,phi)*Ts;

t_conv2 = min(t_conv)+min(t1):Ts:max(t1)+max(t_conv)
%t_conv2 = -2*A*T:Ts:2*A*T+length(Xi_n)-Ts;
f = -Fs/2:Fs/Nf:Fs/2-Fs/Nf; 

plot(t_conv2,Xi_demodulated);
grid on;
title('Output waveform of filtered Xi(t)');
xlabel('t(s)');
ylabel('Xi(t)');

plot(t_conv2,Xq_demodulated);
grid on;
title('Output waveform of filtered Xq(t)');
xlabel('t(s)');
ylabel('Xq(t)');

T_total_i=length(Xi_demodulated)*Ts; %
XF_i=fftshift(fft(Xi_demodulated,Nf))*Ts;  %fourier
PxF_i=(abs(XF_i).^2)/T_total_i %mias ylopoihshs

T_total_q=length(Xq_demodulated)*Ts; %
XF_q=fftshift(fft(Xq_demodulated,Nf))*Ts;  %fourier
PxF_q=(abs(XF_q).^2)/T_total_q %mias ylopoihshs


figure();
plot(f, PxF_i);
title('Periodogram of filtered Xi(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PxF_q);
title('Periodogram of filtered Xq(t)');
xlabel('F');
ylabel('Px(F)');
grid on;

%-----Erwthma 11-----  

Xi_demod_t_sampling = Xi_demodulated((2*A*T/Ts)+1 : over : length(Xi_demodulated)-(2*A*T/Ts));
Xq_demod_t_sampling = Xq_demodulated((2*A*T/Ts)+1 : over : length(Xq_demodulated)-(2*A*T/Ts));

for i=1:N
    Sampling(i,1) = Xi_demod_t_sampling(i);
    Sampling(i,2) = Xq_demod_t_sampling(i);
end

scatterplot(Sampling) 

%-----Erwthma 12-----  

for i=1:N
    Xi_det(i)=detect_4_PAM(Xi_demod_t_sampling(i),A);
    Xq_det(i)=detect_4_PAM(Xq_demod_t_sampling(i),A); 
end

%-----Erwthma 13----- 

%Pinakas sumbolwn pou esteila:
X_sent = [XI ; XQ];

%Pinakas sumbolwn pou antistoixhsa:
X_detect = [Xi_det ; Xq_det];
Matrix= X_sent-X_detect;
error=0;

for i=1:N 
    if(Matrix(1,i)~=0 || (Matrix(2,i)~=0))
    error=error+1;
    end
end

%-----Erwthma 15----- 
Biti=PAM_4_to_bits(Xi_det,A); 
Bitq=PAM_4_to_bits(Xq_det,A);
Bit_est=[Biti Bitq];
Matrix_bit= Bit_est-bit_seq;

errors=0;

for i=1:4*N
    if(Matrix_bit(i) ~=0 )
    errors = errors+1;
    end
end



