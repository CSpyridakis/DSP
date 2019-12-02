% -------------------------------------------------------------------------
%   LAB 3 
%
%   Authors : 
%               - Spyridakis Christos
%               - Marinou Ioanna
%               - Paterakis Isidoros
%
%   Created Date : 27/11/2019
%   Last Updated : 2/12/2019
%
%   Description: 
%               Code created for labs of Digital Signal Processing Course
%               in Technical University of Crete
%
% -------------------------------------------------------------------------

close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1a begins

Fs = 10000;                                    % sampling frequency
Wp = 2*3.14*3000;                              % passband frequency, the cutoff frequency (normalized)
Ws = 2*3.14*5000;                              % stopband frequency (normalized)
Rp = 3;                                        % ripple
Nf = 2048;                                     % samples    
f = linspace(0, Fs/2, Nf);                     % frequency axis
Rs = 30;                                       % attenuation

[n, Wn] = buttord(Wp, Ws, Rp, Rs, 's');        % gets the order and cutoff frequency of the filter
[z, p, k] = buttap(n);                         % gets the zeroes, poles and gain of the butterworth filter
[num, den] = zp2tf(z, p, k);                   % gets the coefficients of the numerator and denominator of the transfer function
[num2, den2] = lp2lp(num, den, Wn);            % changes the frequency from 1 rad/s to Wn and transforms the coefficients of the numerator and denominator
analogfilter1 = freqs(num2, den2, 2*3.14*f);   % analog filter frequency response 
[NUM1, DEN1] = bilinear(num2, den2, Fs);       % converts the the coefficients to the descrete equivalents
digitalfilter1 = freqz(NUM1, DEN1, f, Fs);     % digital filter frequency response

% Display them
figure; hold on;
plot(f,mag2db(abs(analogfilter1)),'--');
plot(f,mag2db(abs(digitalfilter1)));
legend('Analog filter', 'Digital filter'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Comparison between Analog and Digital Butterworth Lowpass Filters.'); ylim([-400 100]); hold off;

%1a ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1b begins

Rs = 50;                                       % change of the attenuation

[n, Wn] = buttord(Wp, Ws, Rp, Rs, 's');        % gets the order and cutoff frequency of the filter
[z, p, k] = buttap(n);                         % gets the zeroes, poles and gain of the butterworth filter
[num, den] = zp2tf(z, p, k);                   % gets the coefficients of the numerator and denominator of the transfer function
[num2, den2] = lp2lp(num, den, Wn);            % changes the frequency from 1 rad/s to Wn and transforms the coefficients of the numerator and denominator
analogfilter2 = freqs(num2, den2, 2*3.14*f);   % analog filter frequency response 
[NUM2, DEN2] = bilinear(num2, den2, Fs);       % converts the the coefficients to the descrete equivalents
digitalfilter2 = freqz(NUM2, DEN2, f, Fs);     % digital filter frequency response

% Display them
figure; hold on;
plot(f,mag2db(abs(analogfilter2)),'--');
plot(f,mag2db(abs(digitalfilter2)));
legend('Analog filter', 'Digital filter'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Comparison between Analog and Digital Butterworth Lowpass Filters.'); ylim([-400 100]); hold off;

%1b ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%2 begins
w = 2;                                      % Wc cutoff frequency
Ts = 0.2;                                   % sampling period
fs = 1/Ts;                                  % sampling frequency
Rp = 3;                                     % ripple
n = [2 16];                                 % order for each case

f=w/(fs*pi);

% Filter
[b1,a1] = cheby1(n(1),Rp,f,'high');
[b2,a2] = cheby1(n(2),Rp,f,'high');

% digital filter frequency response
[h1,w1] = freqz(b1,a1);
[h2,w2] = freqz(b2,a2);

% Display them
figure(); hold on;
plot(w1/pi,20*log10(abs(h1)),'--r')
plot(w2/pi,20*log10(abs(h2)),'-.m')
ylabel('Magnitude [db]'); xlabel('f [* pi rad]'); legend('n=', 'n=16'); title('Chebyshev filter');

%2 ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3a begins

Fs = 10*10^3;                               % Fs=10KHz (Sampling Frequency)
Ts = 1/Fs;

% x(t) sampling
number_of_samples = 500; 
n = 0 : number_of_samples - 1;   
x_a = 1 + cos(1000*n*Ts) + cos(16000*n*Ts) + cos(30000*n*Ts);

% Filters
y1 = filter(NUM1, DEN1, x_a);               % attenuation = 30db
y2 = filter(NUM2, DEN2, x_a);               % attenuation = 50db

% X(F) Fourier Transform of x(t)
F = -Fs/2 : Fs/(number_of_samples) : (Fs/2) - (Fs/number_of_samples);
X_a = abs(fftshift(fft(x_a)*Ts));

% Filters Fourier
Y1 = abs(fftshift(fft(y1)*Ts));
Y2 = abs(fftshift(fft(y2)*Ts));

% Display them
figure(); 
subplot(3,1,1); plot(n, x_a, 'b', n, x_a, 'r.'); legend('Reconstructed','Samples'); legend('Location','NorthEast'); title('x(t)'); xlabel('n'); ylabel('Amplitude');
subplot(3,1,2); plot(n, y1, 'b', n, y1, 'r.'); legend('Reconstructed','Samples'); legend('Location','NorthEast'); title('x(t) filtered attenuation = 30db'); xlabel('n'); ylabel('Amplitude');
subplot(3,1,3); plot(n, y2, 'b', n, y2, 'r.'); legend('Reconstructed','Samples'); legend('Location','NorthEast'); title('x(t) filtered attenuation = 50db'); xlabel('n'); ylabel('Amplitude');

figure(); 
subplot(3,1,1); plot(F, X_a, 'b'); title('Frequency spectrum of x(t) - (|X(F)|)'); xlabel('F'); ylabel('Magnitude');
subplot(3,1,2); plot(F, Y1, 'b'); title('Frequency spectrum of x(t) - (|X(F)|) filtered attenuation = 30db'); xlabel('F'); ylabel('Magnitude');
subplot(3,1,3); plot(F, Y2, 'b'); title('Frequency spectrum of x(t) - (|X(F)|) filtered attenuation = 50db'); xlabel('F'); ylabel('Magnitude');
    
%3a end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3b begins
 
Ts = 0.2;                                   % sampling period
Fs = 1/Ts;                                  % sampling frequency

% x(t) sampling
number_of_samples = 500; 
n = 0 : number_of_samples - 1;   
x_b = 1 + cos(1.5*n*Ts) + cos(5*n*Ts);

% Filter signal
y3 = filter(b2,a2,x_b); 

% X(F) Fourier Transform of x(t)
F = -Fs/2 : Fs/(number_of_samples) : (Fs/2) - (Fs/number_of_samples);
X_b = abs(fftshift(fft(x_b)*Ts));

% Filter's Fourier
Y3 = abs(fftshift(fft(y3)*Ts));

% Display them
figure(); 
subplot(2,1,1); plot(n, x_b, 'b', n, x_b, 'r.'); legend('Reconstructed','Samples'); legend('Location','NorthEast'); title('x(t)'); xlabel('n'); ylabel('Amplitude');
subplot(2,1,2); plot(n, y3, 'b', n, y3, 'r.'); legend('Reconstructed','Samples'); legend('Location','NorthEast'); title('x(t) filtered'); xlabel('n'); ylabel('Amplitude');

figure(); 
subplot(2,1,1); plot(F, X_b, 'b'); title('Frequency spectrum of x(t) - (|X(F)|)'); xlabel('F'); ylabel('Magnitude');
subplot(2,1,2); plot(F, Y3, 'b'); title('Frequency spectrum of x(t) - (|X(F)|) filtered'); xlabel('F'); ylabel('Magnitude');


%3b ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


