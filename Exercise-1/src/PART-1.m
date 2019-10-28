close all
clear all

%1A starts

n = 0:9;       %the time space

signal1 =  zeros(size(n));
signal1(n>4) = 1;              %making the first sequence

figure
stem(n, signal1);               %drawing the first sequence
ylabel('Signal1[n]');                         
xlabel('----->n');
title('The first sequence.');

signal2 = 0;

for k = 1:length(n)
if n(k)<5
signal2(k) = 0;                     %making the second sequence
else
signal2(k) = n(k);
end
end

figure
stem(n, signal2);               %drawing the second sequence
ylabel('Signal2[n]');                         
xlabel('----->n');
title('The second sequence.');

m=length(signal1);                            %saving the length of each sequence
n=length(signal2);

X=[signal1,zeros(1,n)];                     % filling with zeros based on the other sequence's length
H=[signal2,zeros(1,m)];                    

for i=1:n+m-1                           %starts the convolution iteration
    Y(i)=0;
    for t=1:i
        Y(i)=Y(i)+X(t)*H(i-t+1);        %calculates the convolution at the i index for all the different t values
    end
end

figure
stem(Y);                                    %draws the manually made convolution
ylabel('Y[n]');                         
xlabel('----->n');
title('Convolution of two Signals without conv function.');
  
figure
tmp = conv(signal1, signal2);
stem(tmp);                                      %draws the automated convolution
ylabel('Y`[n]');
xlabel('----->n');
title('Convolution of two Signals with conv function.');

%1A ends
%1B starts

signalfur1 = fft(signal1, n+m-1);              %makes the Furrier transform of the 2 sequences
signalfur2 = fft(signal2, n+m-1); 

temp = signalfur1.*signalfur2;          %calculates the product of their Furrier transforms
result = ifft(temp);                    %makes the inverse Furrier of the product

figure
stem(result);                           %draws the results
ylabel('Result[n]');
xlabel('----->n');
title('Convolution result with the use of Furrier Transform.');    
    
    
    
    
    
    
    
    
    
    
    