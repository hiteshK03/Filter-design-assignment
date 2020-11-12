close all;
clearvars;
clc;

% m is unique for everyone
% find addend such that specs are met by looking at plot
% remove semicolon to see the value printed

% 1. Specs
m = 36;
q_m = floor(0.1*m);
r_m = m - 10*q_m;
BL = 25+1.9*q_m + 4.1*r_m;
BH = BL + 20;
trans_bw = 4*10^3;

% 2. Band Edge specifications
fp1 = BL*10^3-trans_bw;
fs1 = BL*10^3;
fs2 = BH*10^3;
fp2 = BH*10^3+trans_bw;
f_samp = 260*10^3;
wp1_by_pi=fp1*2/f_samp;
ws1_by_pi=fs1*2/f_samp;
ws2_by_pi=fs2*2/f_samp;
wp2_by_pi=fp2*2/f_samp;

% 3. Transformed Band Edge specs using Bilinear Transformation         
Wp1 = tan(fp1/f_samp*pi);
Ws1 = tan(fs1/f_samp*pi);
Ws2 = tan(fs2/f_samp*pi);
Wp2 = tan(fp2/f_samp*pi);
B=Wp2-Wp1;
W0=sqrt(Wp1*Wp2);

% 4. Frequency transformation 
WLp1=(B*Wp1)/(W0^2-Wp1^2);
WLs1=(B*Ws1)/(W0^2-Ws1^2);
WLs2=(B*Ws2)/(W0^2-Ws2^2);
WLp2=(B*Wp2)/(W0^2-Wp2^2);

% 5. Lowpass specs
Ws=min(abs(WLs1),abs(WLs2));
Wp=WLp1;

% 6. Chebyshev LPF parameters
delta = 0.15;
D1 = 1/((1-delta)^2)-1;    
D2 = 1/(delta^2)-1;
epsilon = sqrt(D1);
N = ceil(acosh(sqrt(D2)/epsilon)/acosh(Ws/Wp));
% plot poles
sk=zeros(2*N,1);
Bk=asinh(1/epsilon)/N;
for k=1:(2*N)
    sk(k)=-Wp*sin(pi*(2*k-1)/(2*N))*sinh(Bk)+j*Wp*cos(pi*(2*k-1)/(2*N))*cosh(Bk);
end
figure(1),plot(real(sk),imag(sk),'rX')  % poles
title("Poles of the analog lowpass transfer function")
xlabel("Re(s)")
ylabel("Im(s)")
axis equal
grid on
sk_maj = Wp*cosh(Bk);
sk_min=Wp*sinh(Bk);
t = linspace(0,2*pi,100);
hold on
plot(sk_min*cos(t),sk_maj*sin(t),'b-')  % circle of radius Wc

% 8. bandpass transfer function
numerator=sk(1)*sk(2)*sk(3)*sk(4)/sqrt(1+epsilon^2);
[num,den] = zp2tf([],[sk(1:4)],numerator);   %TF with poles sk and numerator Wc^N and no zeroes
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(W0*W0+s*s));        %bandstop transformation
%coeffs of analog bpf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                         %collect coeffs into matrix form
kd= ds(1); 
kn = ns(1);
ds = ds/kd;
ns = ns/kn;
k=kn/kd;

% 9. coeffs of discrete BPF
discrete_bsf(z) = analog_bsf((1-z)/(z+1));              %bilinear transformation
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
kd = dz(1); %normalisation factor
kn = nz(1); %normalisation factor
dz = dz/kd;
nz = nz/kn;
k=kn/kd;
nz=nz*k;
fvtool(nz,dz,'Analysis','freq');                                       %frequency response in dB

[~,p,~]=tf2zp(nz,dz);
figure;
plot(real(p),imag(p),'rX');
title("Poles of the transfer function")
xlabel("Re(z)")
ylabel("Im(z)")
axis equal
grid on
t = linspace(0,2*pi,1000);
hold on
plot(cos(t),sin(t),'b-') 

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")
xlabel("Hz")
ylabel("|H(f)|")
xline(fs1,'--m');
xline(fp1,'--m');
xline(fp2,'--m');
xline(fs2,'--m');
yline(0.85,'--m');
yline(0.15,'--m');
grid