%Ex 2 begins %%%%%%%%%%

close all
clear all

numerator = [4 -3.5 0];           %adds the factors of the numerator to a vector
denominator = [1 -2.5 1];           %adds the factors of the denominator to a vector
[roots, poles, k] = residuez(numerator, denominator)   %uses the residuez MATLAB function to calculate the roots, poles and k(not important here)

syms z;                 %the definition of z
        
h1 = roots(1) / ((1-poles(1)*z^(-1)));      %makes the first half of the analyzed transfer function
h2 = roots(2) / ((1-poles(2)*z^(-1)));      %makes the second half of the analyzed transfer function
H = h1 + h2;                                %adds up the 2 parts
pretty(H)                                   %prints the result

Hz = iztrans(H)                               %calculates the inverse z transform of the transfer function