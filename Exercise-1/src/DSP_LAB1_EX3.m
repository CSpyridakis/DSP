%A
close all
clear all

N=128;
T=0.05;
tmax = T;


t= linspace(-tmax,tmax,N);
x = 10*cos(2*pi*20*t)-4*sin(2*pi*40*t+5);		%Normal signal

figure;
plot(t,x)										%Designing signal
hold on;
Ts = 2*T/N;										%Period of sampling
td = -tmax:Ts:tmax;								
x1 = 10*cos(2*pi*20*td)-4*sin(2*pi*40*td+5);	%Sampled signal
plot(td,x1,'r.');								%Designing sampled signal (with read dots)
%Adding details to the diagram
title('Signal and samples');
xlabel('t(s)');
ylabel('x(t)');
hold off;

figure;										%Creating new figure
Fs=1/Ts;										%
X = fftshift(fft(x1,N)*Ts);						%Creating X(F) from the sampled signal
F = [-Fs/2:Fs/N:Fs/2-Fs/N];						%F vector declaration
plot(F,abs(X));									%Plotting the Fourier transformation
%Adding details to the diagram
xlim([-100,100]);
title('Spectrum of the given function');
xlabel('F(Hz)');
ylabel('X(F)');


%B
T=0.05;
N = 2048;
tmax = T;
l = randi(30);									%Generating random int to add to the signal
Fs = 8000;
Ts = 1/Fs;
td = -tmax:Ts:tmax;

Fs=1/Ts;
df = 125;
f = [100:df:475];
for k = 1:length(f)
    figure;									%Creating new figure
    x = sin(2*pi*f(k).*td + l);					%Creating x(t) signal
    X = fftshift(fft(x,N)*Ts);					%Generating Fourier transformation of x(t)
    F = [-Fs/2:Fs/N:Fs/2-Fs/N];					%F vector declaration
    plot(F,abs(X));								%Plotting X(F)
	%Adding details to the diagram
    title(['Spectrum of signal when f = ' num2str(f(k))]);
    xlabel('F(Hz)');
    ylabel('X(F)');
    xlim([-600 600]);
end
T=1/8000;
td = -tmax:Ts:tmax;
f1 = [7500:df:7900];
for k = 1:length(f1)
    figure;                             		%Creating new figure
    x1 = sin(2*pi*f1(k).*td + l);				%Creating x(t) signal
    X1 = fftshift(fft(x1,N)*Ts);				%Generating Fourier transformation
    F = [-Fs/2:Fs/N:Fs/2-Fs/N];					%F vector declaration
    plot(F,abs(X1));							%Plotting X(F)
	%Adding details to the diagram
    title(['Spectrum of signal when f = ' num2str(f1(k))]); 
    xlabel('F(Hz)');
    ylabel('X(F)');
    xlim([-600 600]);
end