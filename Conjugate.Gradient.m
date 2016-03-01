%{ 
Jacob Leonard
MATH 467 - Fall 2015
jaleonar@usc.edu
Revision History
Date             Changes           Programmer 
-------------------------------------------------
11/6/2015        Original          Jacob Leonard

%}

%this script is for the conjugate gradient method with Fletcher-Reeves formula

%determine x(0) and y(0) for the start of the methods. The script will
%eventually loop through all of these points in the 100 by 100 square of
%values, and evaluate the number of steps, as well as the rate of
%convergence
for j = 1:100
    x(j) = (-2)+((4*j)/100);
    y(j) = (-2)+((4*j)/100);
end

syms x y 
f = (x^4 + y^4 -6*x^2*y^2-1)^2+(4*x^3*y-4*x*y^3)^2;
X = diff(f,x);
Y = diff(f,y);
XX = diff(X,x);
YY = diff(Y,y);
XY = diff(
YX = diff(
YX = diff(
