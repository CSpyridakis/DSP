%-------------------------------------------------------------
%
%
%
%
%-------------------------------------------------------------
clear all
close all
clc

% Just for saving in a separate folder figures as png files
DEBUG = true; dirpath = 'photos'; if ~DEBUG && ~exist(dirpath,'dir') ; mkdir(dirpath); end

%
%   Exercise 2
%

% 2.1
% Initial signal 
tmin = 0 ; tmax = 0.5 ;
dt = 0.001;                             % dt = 1ms
t = [ tmin + dt : dt : tmax - dt ];     % 0 < t < 500ms
x = 5*cos(24*pi*t)-2*sin(1.5*pi*t);     % x(t) = 5cos(24πt)-2sin(1.5πt)
% Plot Initial signal
stepName = '2.1 Initial signal';
f = figure();
p1 = plot(t,x,'b'); 
legend([p1],'x(t)'); legend('Location','NorthEast');
title(stepName); ylabel('x(t)'); xlabel('t(sec)'); %axis([xmin xmax ymin ymax]);
if ~DEBUG ; saveas(f,strcat(dirpath, '/', stepName, '.jpg')) ; end

% 2.2
% Fourier Tranform
NFFT = 2^nextpow2(length(x));       % NFFT = biggest power of 2 
X=fftshift(fft(x,NFFT));

% 2.3
% Sampling Frequences
Ts = [1/48 1/24 1/12];
Ts_str = ['48';'24';'12'];

stepName = '2.3 ';
for i=1:length(Ts)
    n = [ceil(tmin/Ts(i)) : floor(tmax/Ts(i))];                       % length(n) = nunber of samples
    x_sampled = 5*cos(24*pi*n*Ts(i))-2*sin(1.5*pi*n*Ts(i));
    
    % New figure
    f = figure(); grid on;
    
    % Sampled signal
    subplot(2,1,1)
    p1 = plot(t,x,'b'); 
    hold on;
    p2 = plot(n*Ts(i),x_sampled,'r.');
    legend([p1,p2],'x(t)','X(s)'); legend('Location','NorthEast');
    title(strcat(stepName, 'Sampled - Ts=1/', Ts_str(i,:))); ylabel('x(t)'); xlabel('t(sec)'); %axis([0 .5 4 -4]);
    
    % Reconstructed signal
    subplot(2,1,2)
    p1 = plot(t,x,'b'); 
    hold on;
    p2 = plot(n*Ts(i),x_sampled,'r');
    legend([p1,p2],'x(t)','X(s)'); legend('Location','NorthEast');
    title(strcat(stepName, 'Reconstructed - Ts=1/', Ts_str(i,:))); ylabel('x(t)'); xlabel('t(sec)'); %axis([0 .5 4 -4]);
    
    if ~DEBUG ; saveas(f,strcat(dirpath, '/', stepName, '- Ts is 1_', Ts_str(i,:), '.jpg')) ; end 
end








