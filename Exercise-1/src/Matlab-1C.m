%A
N=128;
T=0.05;
tmax = T;

t= linspace(-tmax,tmax,N);
x = 10*cos(2*pi*20*t)+4*sin(2*pi*40*t+5);
plot(t,x)
hold on;
Ts = 2*T/N;
td = -tmax:Ts:tmax;
x1 = 10*cos(2*pi*20*td)+4*sin(2*pi*40*td+5);
plot(td,x1,'r.');
title("Signal and samples");
xlabel("t(s)");
ylabel("x(t)");


figure(2);
Fs=1/Ts;
X = fftshift(fft(x1,N)*Ts);
F = [-Fs/2:Fs/N:Fs/2-Fs/N];
plot(F,abs(X));
xlim([-100,100]);
title("Spectrum of the given function");
xlabel("F(Hz)");
ylabel("X(F)");


%B
T=0.05;
N = 2048;
tmax = T;
l = randi(30);
Fs = 8000;
Ts = 1/Fs;
td = -tmax:Ts:tmax;

Fs=1/Ts;
df = 125;
f = [100:df:475];
for k = 1:length(f)
    figure(k);
    x = sin(2*pi*f(k).*td + l);
    X = fftshift(fft(x,N)*Ts);
    F = [-Fs/2:Fs/N:Fs/2-Fs/N];
    plot(F,abs(X));
    title(["Spectrum of signal when f = "+f(k)]);
    xlabel("F(Hz)");
    ylabel("X(F)");
    xlim([-600 600]);
end
T=1/8000;
td = -tmax:Ts:tmax;
f1 = [7500:df:7900];
for k = 1:length(f1)
    figure(k+length(f));
    x1 = sin(2*pi*f1(k).*td + l);
    X1 = fftshift(fft(x1,N)*Ts);
    F = [-Fs/2:Fs/N:Fs/2-Fs/N];
    plot(F,abs(X1));
    title("Spectrum of signal when f = "+[f1(k)]); 
    xlabel("F(Hz)");
    ylabel("X(F)");
    xlim([-600 600]);
end