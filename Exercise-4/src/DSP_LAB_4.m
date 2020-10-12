% -------------------------------------------------------------------------
%   LAB 4 
%
%   Authors : 
%               - Spyridakis Christos
%               - Marinou Ioanna
%               - Paterakis Isidoros
%
%   Created Date : 6/11/2019
%   Last Updated : 16/12/2019
%
%   Description: 
%               Code created for labs of Digital Signal Processing Course
%               in Technical University of Crete
%
% -------------------------------------------------------------------------

close all; clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1 begins
Wc = 0.4 * pi;
Fs = 100;               % Sampling Frequency
N = 21;                 % Window Length

Fc = Wc / 2 * pi;
Wn = Fc / (Fs/2);       % Normalized Cutoff Frequency

% Create Rectangular and Hamming windows
rect_window = rectwin(N);
hamm_window = hamming(N);

% Create FIR filters using the window method
rect_filter = fir1(N-1, Wn, rect_window);
hamm_filter = fir1(N-1, Wn, hamm_window);

% Find out frequency response
[Hrect, Wrect] = freqz(rect_filter, N);
[Hhamm, Whamm] = freqz(hamm_filter, N);

% Plot them
figure()
plot(Wrect, abs(Hrect), 'b', Whamm, abs(Hhamm), 'r'); legend('Rectangular','Hamming');
title('Rectangular and Hamming Frequency Responses'); xlabel('Frequency'); ylabel('Amplitute'); 

%1 ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%B1 begins

Wc = 0.5 * pi;              % cutoff frequency
Fc = Wc / (2*pi)
Fs = 100;                   % sampling frequency
N1 = 21;                    % first window length
N2 = 41;                    % second window length
Wn = Fc / (Fs/2);           % normalized cutoff frequency

hammingn1 = hamming(N1);    % makes the first hamming window
hammingn2 = hamming(N2);    % makes the second hamming window

hanningn1 = hann(N1);       % makes the first hanning window
hanningn2 = hann(N2);       % makes the second hanning window

hammingfilter1 = fir1(N1-1, Wn, hammingn1);     % makes the first window based hamming FIR digital filter
hammingfilter2 = fir1(N2-1, Wn, hammingn2);     % makes the second window based hamming FIR digital filter
hanningfilter1 = fir1(N1-1, Wn, hanningn1);     % makes the first window based hanning FIR digital filter
hanningfilter2 = fir1(N2-1, Wn, hanningn2);     % makes the second window based hanning FIR digital filter

[hhamm1,whamm1] = freqz(hammingfilter1, N1);    % calculates the frequency response of the first hamming filter
[hhamm2,whamm2] = freqz(hammingfilter2, N2);    % calculates the frequency response of the second hamming filter
[hhann1,whann1] = freqz(hanningfilter1, N1);    % calculates the frequency response of the first hanning filter
[hhann2,whann2] = freqz(hanningfilter2, N2);    % calculates the frequency response of the second hanning filter

HHAMM1 = abs(hhamm1);   % calculates the absolute value of the first hamming filter frequency response
HHAMM2 = abs(hhamm2);   % calculates the absolute value of the second hamming filter frequency response
HHANN1 = abs(hhann1);   % calculates the absolute value of the first hanning filter frequency response
HHANN2 = abs(hhann2);   % calculates the absolute value of the second hanning filter frequency response

figure;
subplot(1,2,1);
plot(whamm1, HHAMM1);
xlabel('Frequency'); ylabel('Magnitude'); title('Hamming: N = 21');
subplot(1,2,2);                 %prints the results of the first half
plot(whamm2, HHAMM2); 
xlabel('Frequency'); ylabel('Magnitude'); title('Hamming: N = 41');

figure;
subplot(1,2,1);
plot(whann1, HHANN1);
xlabel('Frequency'); ylabel('Magnitude'); title('Hanning: N = 21');
subplot(1,2,2);                 %prints the results of the second half
plot(whann1, HHANN2);
xlabel('Frequency'); ylabel('Magnitude'); title('Hanning: N = 41');

%B1 ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%B2 begins
N = 500;
Tmax = 2;
Fs = 100;
Ts = 1/Fs;
t = [0:N-1]*Ts;
f = [-Fs/2:Fs/N:Fs/2-Fs/N];
x = sin(15*t) + 0.25*sin(200*t);

%Filtering the signal through the four different filters
y1 = filter(hammingfilter1, 1, x);  
y2 = filter(hammingfilter2, 1, x);
y3 = filter(hanningfilter1, 1, x);
y4 = filter(hanningfilter2, 1, x);

%Creating the Fourier transformation for the for different signals
Y1 = abs(fftshift(fft(y1)*Ts));
Y2 = abs(fftshift(fft(y2)*Ts));
Y3 = abs(fftshift(fft(y3)*Ts));
Y4 = abs(fftshift(fft(y4)*Ts));
X = abs(fftshift(fft(x,N)*Ts));

%Plotting the signals
figure();
subplot(2,1,1);
plot(f,X); hold on;
plot(f,Y1);
title(strcat('Fs = 100 Hamming filter N = ', num2str(N1))); legend('X(f)','Y(f)');

subplot(2,1,2);
plot(f,X); hold on;
plot(f,Y2);
title(strcat('Fs = 100 Hamming filter N = ', num2str(N2))); legend('X(f)','Y(f)');

figure();
subplot(2,1,1);
plot(f,X); hold on;
plot(f,Y3);
title(strcat('Fs = 100 Hanning filter N = ', num2str(N1))); legend('X(f)','Y(f)');

subplot(2,1,2);
plot(f,X); hold on;
plot(f,Y4);
title(strcat('Fs = 100 Hanning filter N = ', num2str(N2))); legend('X(f)','Y(f)');
%B2 ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%B3 begins

Fs = 50;
Ts = 1/Fs;
t = [0:N-1]*Ts;
x = sin(15*t) + 0.25*sin(200*t);
f = [-Fs/2:Fs/N:Fs/2-Fs/N];

%Filtering the signal through the four different filters
y1 = filter(hammingfilter1, 1, x);
y2 = filter(hammingfilter2, 1, x);
y3 = filter(hanningfilter1, 1, x);
y4 = filter(hanningfilter2, 1, x);

%Creating the Fourier transformation for the for different signals
Y1 = abs(fftshift(fft(y1)*Ts));
Y2 = abs(fftshift(fft(y2)*Ts));
Y3 = abs(fftshift(fft(y3)*Ts));
Y4 = abs(fftshift(fft(y4)*Ts));
X = abs(fftshift(fft(x,N)*Ts));

%Plotting the signals
figure();
subplot(2,1,1);
plot(f,X); hold on;
plot(f,Y1);
title(strcat('Fs = 50 Hamming filter N = ', num2str(N1))); legend('X(f)','Y(f)');

subplot(2,1,2);
plot(f,X); hold on;
plot(f,Y2);
title(strcat('Fs = 50 Hamming filter N = ', num2str(N2))); legend('X(f)','Y(f)');

figure();
subplot(2,1,1);
plot(f,X); hold on;
plot(f,Y3);
title(strcat('Fs = 50 Hanning filter N = ', num2str(N1))); legend('X(f)','Y(f)');

subplot(2,1,2);
plot(f,X); hold on;
plot(f,Y4);
title(strcat('Fs = 50 Hanning filter N = ', num2str(N2))); legend('X(f)','Y(f)');

%B3 ends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
