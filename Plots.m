%{
Jacob Leonard
MATH 467 - Fall 2015
jaleonar@usc.edu
Revision History
Date             Changes           Programmer
-------------------------------------------------
11/6/2015        Original          Jacob Leonard

%}

%this script is for general graphs and charts, and evaluating the function
%at all the given points

%define the values of x and y from -2 to 2, increasing by 1/25, for 100
%values
for j = 1:101
    x(j) = (-2)+((4*(j-1))/100);
    y(j) = (-2)+((4*(j-1))/100);
end

%these plots and functions are for the complex function

%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(x,y) (abs((x+1i*y)^4-1)^2);
G = {@(x,y) (8*((x+1i*y)^4-1)*(x+1i*y)^2);@(x,y) (8*1i*((x+1i*y)^4-1)*(x+1i*y)^2)};
%Gradient = [g{1}(x,y),g{2}(x,y)];
%when the
H = {@(x,y) (24*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*1i*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*1i*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*(x+1i*y)^2*(-(7/3)*(x+1i*y)^4+1))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%plot the original function on its own
for i = 1:101
    for j = 1:101
        Values(i,j) = f(x(i),y(j));
    end
end
%{
contourf(x,y,Values);
xlabel('X');
ylabel('Y');
title('Function Values Evaluated at X=[-2:2], Y=[-2:2]');

%plot the expanded function real values
for i = 1:101
    for j = 1:101
        Values1(i,j) = f(x(i),y(j));
    end
end
%plot the real values
V1 = real(Values1);
contourf(x,y,V1);
xlabel('X');
ylabel('Y');
%zlabel('Real Values f(X,Y)');
title('Real Function Values Evaluated at X=[-2:2], Y=[-2:2]');

%plot the imaginary values
V2 = imag(Values1);
contourf(x,y,V2);
xlabel('X');
ylabel('Y');
%zlabel('Imaginary Values f(X,Y)');
title('Imaginary Function Values Evaluated at X=[-2:2], Y=[-2:2]');
%}
%compute the gradient of the function along the roots
for i = 1:101
    G1(:,:,1) = [G{1}(x(i),1i*x(i)+1),G{2}(x(i),1i*x(i)+1)];
    G2(:,:,1) = [G{1}(x(i),1i*x(i)-1),G{2}(x(i),1i*x(i)-1)];
    G4(:,:,1) = [G{1}(x(i),1i*x(i)-1i),G{2}(x(i),1i*x(i)-1i)];
    G3(:,:,1) = [G{1}(x(i),1i*x(i)+1i),G{2}(x(i),1i*x(i)+1i)];
end

%compute the hessian of the function along the roots
for i = 1:101
    H1(:,:,i) = [H{1}(x(i),1i*x(i)+1),H{2}(x(i),1i*x(i)+1);H{3}(x(i),1i*x(i)+1),H{4}(x(i),1i*x(i)+1)];
    H2(:,:,i) = [H{1}(x(i),1i*x(i)-1),H{2}(x(i),1i*x(i)-1);H{3}(x(i),1i*x(i)-1),H{4}(x(i),1i*x(i)-1)];
    H4(:,:,i) = [H{1}(x(i),1i*x(i)-1i),H{2}(x(i),1i*x(i)-1i);H{3}(x(i),1i*x(i)-1i),H{4}(x(i),1i*x(i)-1i)];
    H3(:,:,i) = [H{1}(x(i),1i*x(i)+1i),H{2}(x(i),1i*x(i)+1i);H{3}(x(i),1i*x(i)+1i),H{4}(x(i),1i*x(i)+1i)];
end

%compute the gradient of the function near the roots, delta = .04
for i = 1:101
    G1addD(:,:,i) = [G{1}(x(i),1i*x(i)+1.04),G{2}(x(i),1i*x(i)+1.04)];
    G2addD(:,:,i) = [G{1}(x(i),1i*x(i)-.96),G{2}(x(i),1i*x(i)-.96)];
    G4addD(:,:,i) = [G{1}(x(i),1i*x(i)-1i+.04),G{2}(x(i),1i*x(i)-1i+.04)];
    G3addD(:,:,i) = [G{1}(x(i),1i*x(i)+1i+.04),G{2}(x(i),1i*x(i)+1i+.04)];
end

%compute the gradient of the function near the roots, delta = -.04
for i = 1:101
    G1minusD(:,:,i) = [G{1}(x(i),1i*x(i)+.96),G{2}(x(i),1i*x(i)+.96)];
    G2minusD(:,:,i) = [G{1}(x(i),1i*x(i)-1.04),G{2}(x(i),1i*x(i)-1.04)];
    G4minusD(:,:,i) = [G{1}(x(i),1i*x(i)-1i-.04),G{2}(x(i),1i*x(i)-1i-.04)];
    G3minusD(:,:,i) = [G{1}(x(i),1i*x(i)+1i-.04),G{2}(x(i),1i*x(i)+1i-.04)];
end

%compute the hessian of the function along the roots
for i = 1:101
    Hadd1(:,:,i) = [H{1}(x(i),1i*x(i)+1.04),H{2}(x(i),1i*x(i)+1.04);H{3}(x(i),1i*x(i)+1.04),H{4}(x(i),1i*x(i)+1.04)];
    Hadd2(:,:,i) = [H{1}(x(i),1i*x(i)-.96),H{2}(x(i),1i*x(i)-.96);H{3}(x(i),1i*x(i)-.96),H{4}(x(i),1i*x(i)-.96)];
    Hadd4(:,:,i) = [H{1}(x(i),1i*x(i)-1i+.04),H{2}(x(i),1i*x(i)-1i+.04);H{3}(x(i),1i*x(i)-1i+.04),H{4}(x(i),1i*x(i)-1i+.04)];
    Hadd3(:,:,i) = [H{1}(x(i),1i*x(i)+1i+.04),H{2}(x(i),1i*x(i)+1i+.04);H{3}(x(i),1i*x(i)+1i+.04),H{4}(x(i),1i*x(i)+1i+.04)];
end

%compute the hessian of the function along the roots
for i = 1:101
    Hminus1(:,:,i) = [H{1}(x(i),1i*x(i)+.96),H{2}(x(i),1i*x(i)+.96);H{3}(x(i),1i*x(i)+.96),H{4}(x(i),1i*x(i)+.96)];
    Hminus2(:,:,i) = [H{1}(x(i),1i*x(i)-1.04),H{2}(x(i),1i*x(i)-1.04);H{3}(x(i),1i*x(i)-1.04),H{4}(x(i),1i*x(i)-1.04)];
    Hminus4(:,:,i) = [H{1}(x(i),1i*x(i)-1i-.04),H{2}(x(i),1i*x(i)-1i-.04);H{3}(x(i),1i*x(i)-1i-.04),H{4}(x(i),1i*x(i)-1i-.04)];
    Hminus3(:,:,i) = [H{1}(x(i),1i*x(i)+1i-.04),H{2}(x(i),1i*x(i)+1i-.04);H{3}(x(i),1i*x(i)+1i-.04),H{4}(x(i),1i*x(i)+1i-.04)];
end

%these plots and charts are for the real valued function

%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(x,y) ((x^4+y^4-6*x^2*y^2-1)^2+(4*x^3*y-4*x*y^3)^2);
G = {@(x,y) (8*x*(x^6+3*x^4*y^2+x^2*(3*y^4-1)+y^2*(y^4+3))),@(x,y) (8*y*(x^6+3*x^4*y^2+3*x^2*(y^4+1)+y^2*(y^4-1)))};
%Gradient = [g{1}(x,y),g{2}(x,y)];
%when the
H = {@(x,y) 8*(7*x^6+15*x^4*y^2+x^2*(9*y^4-3)+y^2*(y^4+3)),@(x,y) 48*x*y*(x^4+2*x^2*y^2+y^4+1);@(x,y) 48*x*y*(x^4+2*x^2*y^2+y^4+1),@(x,y) 8*(x^6+9*x^4*y^2+3*x^2*(5*y^4+1)+y^2*(7*y^4-3))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%define the values of x and y from -2 to 2, increasing by 1/25, for 100
%values
for j = 1:101
    x(j) = (-2)+((4*(j-1))/100);
    y(j) = (-2)+((4*(j-1))/100);
end


%plot the original function on its own
for i = 1:101
    for j = 1:101
        Values(i,j) = f(x(i),y(j));
    end
end
%{
contourf(x,y,Values);
xlabel('X');
ylabel('Y');
title('Function Values Evaluated at X=[-2:2], Y=[-2:2]');
%}

%compute the gradient of the function along the roots
G1(:,:) = [G{1}(1,0),G{2}(1,0)];
G2(:,:) = [G{1}(-1,0),G{2}(-1,0)];
G3(:,:) = [G{1}(0,1),G{2}(0,1)];
G4(:,:) = [G{1}(0,-1),G{2}(0,-1)];


%compute the hessian of the function along the roots
H1(:,:) = [H{1}(1,0),H{2}(1,0);H{3}(1,0),H{4}(1,0)];
H2(:,:) = [H{1}(-1,0),H{2}(-1,0);H{3}(-1,0),H{4}(-1,0)];
H3(:,:) = [H{1}(0,1),H{2}(0,1);H{3}(0,1),H{4}(0,1)];
H4(:,:) = [H{1}(0,-1),H{2}(0,-1);H{3}(0,-1),H{4}(0,-1)];


%compute the gradient of the function near the roots, delta = .04
G1addD(:,:) = [G{1}(1,0.04),G{2}(1,0.04)];
G2addD(:,:) = [G{1}(-1,0.04),G{2}(-1,0.04)];
G3addD(:,:) = [G{1}(0,1.04),G{2}(0,1.04)];
G4addD(:,:) = [G{1}(0,-.96),G{2}(0,-.96)];


%compute the gradient of the function near the roots, delta = -.04
G1minusD(:,:) = [G{1}(1,-0.04),G{2}(1,-0.04)];
G2minusD(:,:) = [G{1}(-1,-0.04),G{2}(-1,-0.04)];
G3minusD(:,:) = [G{1}(0,.96),G{2}(0,.96)];
G4minusD(:,:) = [G{1}(0,-1.04),G{2}(0,-1.04)];


%compute the hessian of the function along the roots, delta= +.04
Hadd1(:,:) = [H{1}(1,0.04),H{2}(1,0.04);H{3}(1,0.04),H{4}(1,0.04)];
Hadd2(:,:) = [H{1}(-1,0.04),H{2}(-1,0.04);H{3}(-1,0.04),H{4}(-1,0.04)];
Hadd3(:,:) = [H{1}(0,1.04),H{2}(0,1.04);H{3}(0,1.04),H{4}(0,1.04)];
Hadd4(:,:) = [H{1}(0,-.96),H{2}(0,-.96);H{3}(0,-.96),H{4}(0,-.96)];


%compute the hessian of the function along the roots, delta= -.04

Hminus1(:,:) = [H{1}(1,-0.04),H{2}(1,-0.04);H{3}(1,-0.04),H{4}(1,-0.04)];
Hminus2(:,:) = [H{1}(-1,-0.04),H{2}(-1,-0.04);H{3}(-1,-0.04),H{4}(-1,-0.04)];
Hminus3(:,:) = [H{1}(0,.96),H{2}(0,.96);H{3}(0,.96),H{4}(0,.96)];
Hminus4(:,:) = [H{1}(0,-1.04),H{2}(0,-1.04);H{3}(0,-1.04),H{4}(0,-1.04)];


%this section is for the updated z funciton
%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(z1,z2) (z1^2+4*z1^(3/2)*z2^(1/2)+6*z1*z2-2*z1+4*z1^(1/2)*z2^(3/2)+12*z1^(1/2)*z2^(1/2)+z2^2-2*z2+1);
G = {@(z1,z2) (2*z1+6*z1^(1/2)*z2^(1/2)+6*z2-2+2*z1^(-1/2)*z2^(3/2)+6*z1^(-1/2)*z2^(1/2)),@(x,y) (2*z2+6*z2^(1/2)*z1^(1/2)+6*z1-2+2*z2^(-1/2)*z1^(3/2)+6*z2^(-1/2)*z1^(1/2))};
%Gradient = [g{1}(x,y),g{2}(x,y)];
%when the
H = {@(x,y) (2+3*z1^(-1/2)*z2^(1/2)-z1^(-3/2)*z2^(3/2)-3*z1^(-3/2)*z2^(1/2)),@(x,y) (3*z1^(1/2)*z2^(-1/2)+6+3*z1^(-1/2)*z2^(1/2)+3*z1^(-1/2)*z2^(-1/2));@(x,y) (3*z1^(1/2)*z2^(-1/2)+6+3*z1^(-1/2)*z2^(1/2)+3*z1^(-1/2)*z2^(-1/2)),@(x,y) (2+3*z2^(-1/2)*z1^(1/2)-z2^(-3/2)*z1^(3/2)-3*z2^(-3/2)*z1^(1/2))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%these plots are for the use of z1 and z2 to develop the quadratic function
for j = 1:101
    x(j) = ((2*(j-1))/100);
    y(j) = ((2*(j-1))/100);
end

z1 = x.^4;
z2 = y.^4;

for i = 1:101
    for j = 1:101
        ZValues(i,j) = f(z1(i),z2(j));
    end
end

%this is a general plot for the z function 
contourf(z1,z2,ZValues)
xlabel('z1');
ylabel('z2');
title('Function Values Evaluated at z1=[0:16], z2=[0:16]');

