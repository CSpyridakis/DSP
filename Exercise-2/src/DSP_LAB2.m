% ---------------------------------------------------------------------------------
%   LAB 2 
%
%   Authors : 
%               - Spyridakis Christos
%               - Marinou Ioanna
%               - Paterakis Isidoros
%
%   Created Date : 2/11/2019
%   Last Updated : 5/11/2019
%
%   Description: 
%               Code created for labs of Digital Signal Processing Course
%               in Technical University of Crete
%
% ---------------------------------------------------------------------------------
close all ; clear all ; clc ; 

%---------------------------------------------------------
%    Ex 1 begins 
%-------------------
Ts=1;                           % Sampling frequency

% b
num_1b = [0 0.2 0];             % numetator of transfer function 
den_1b = [1 -0.7 -0.18];        % denominator of transfer function
system_1b=tf(num_1b,den_1b);    % Convertion to Transfer function

figure()
zplane(num_1b,den_1b);          % Display Zero-Pole map including unit circle
title('Zero-Pole Map');

% Second option using system, but without unit cicle display
% figure()
% pzmap(system_1b)


%d

%e


%---------------------------------------------------------
%    Ex 2 begins 
%-------------------
numerator_2a = [4 -3.5 0];                                  %adds the factors of the numerator to a vector
denominator_2a = [1 -2.5 1];                                %adds the factors of the denominator to a vector
[roots, poles, k] = residuez(numerator_2a, denominator_2a)  %uses the residuez MATLAB function to calculate the roots, poles and k(not important here)

syms z;                                                     %the definition of z
        
h1 = roots(1) / ((1-poles(1)*z^(-1)));                      %makes the first half of the analyzed transfer function
h2 = roots(2) / ((1-poles(2)*z^(-1)));                      %makes the second half of the analyzed transfer function
H = h1 + h2;                                                %adds up the 2 parts
pretty(H)                                                   %prints the result

Hz = iztrans(H)                                             %calculates the inverse z transform of the transfer function