% -------------------------------------------------------------------------
%   LAB 2 
%
%   Authors : 
%               - Spyridakis Christos
%               - Marinou Ioanna
%               - Paterakis Isidoros
%
%   Created Date : 2/11/2019
%   Last Updated : 17/11/2019
%
%   Description: 
%               Code created for labs of Digital Signal Processing Course
%               in Technical University of Crete
%
% -------------------------------------------------------------------------
close all ; clear all ; clc ; 

%---------------------------------------------------------
%    Ex 1 begins 
%-------------------
Ts=1;                           % Sampling frequency

% b
num_1b = [0 0.2 0];             % numetator (num) of transfer function 
den_1b = [1 -0.7 -0.18];        % denominator (den) of transfer function
system_1b=tf(num_1b,den_1b);    % Convertion to Transfer function

figure()
zplane(num_1b,den_1b);          % Display Zero-Pole map including unit circle
title('Zero-Pole Map');

% Second option using system, but without unit cicle display
% figure()
% pzmap(system_1b)

%d
df = pi/128;
f = [-pi : df : pi];            
b = [0 0.2 0];                  % num of transfer function  
a = [1 -0.7 -0.18];             % den of transfer function

figure() 
freqz(b,a,f)
title('Diagrams of freqz with n parameter');

figure();
freqz(b,a)
title('Diagrams of freqz without n parameter');

%e
b = [0 0 0.2 0];                % new num of tf with pole added on 1 
a = [1 -1.7 0.52 0.18];         % new den of tf with pole added on 1 

figure(); 
freqz(b,a,f)
title('Diagrams of freqz with an added pole');

%---------------------------------------------------------
%    Ex 2 begins 
%-------------------

% Adds the factors of the numerator to a vector
num_2a = [4 -3.5 0];                                  
% Adds the factors of the denominator to a vector
den_2a = [1 -2.5 1];                                

% Uses the residuez MATLAB function to calculate the roots, poles and 
% k(not important here)
[zero_roots, pole_roots, k] = residuez(num_2a, den_2a)  

syms z;                        %the definition of z
       
%makes the first half of the analyzed transfer function
h1 = zero_roots(1) / ((1-pole_roots(1)*z^(-1)));            
%makes the second half of the analyzed transfer function
h2 = zero_roots(2) / ((1-pole_roots(2)*z^(-1)));            

H = h1 + h2;                    %adds up the 2 parts
pretty(H)                       %prints the result

Hz = iztrans(H) % calculates the inverse z transform of the transfer function