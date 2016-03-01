%{
Jacob Leonard
MATH 467 - Fall 2015
jaleonar@usc.edu
Revision History
Date             Changes           Programmer
-------------------------------------------------
11/6/2015        Original          Jacob Leonard

%}

%this script is for the fixed step size gradient method

%determine x(0) and y(0) for the start of the methods
for i = 1:100
    x(i) = (-2)+((4*i)/100);
    y(i) = (-2)+((4*i)/100);
end

%find the x and y derivatives of the function, and the hessian
%f = ((x^4 + y^4 -6*x^2*y^2-1)^2+(4*x^3*y-4*x*y^3)^2);


%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(x,y) (abs(x^4 + y^4 -6*x^2*y^2-1)^2+1i*(4*x^3*y-4*x*y^3)-1)^2;
G = {@(x,y) (8*((x+1i*y)^4-1)*(x+1i*y)^2);@(x,y) (8*1i*((x+1i*y)^4-1)*(x+1i*y)^2)};
%Gradient = [g{1}(x,y),g{2}(x,y)];
%when the 
H = {@(x,y) (24*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*1i*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*(x+1i*y)^2*(-(7/3)*(x+1i*y)^4+1))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%desired level of accuracy
tolerance = 10^(-7);

%this matrix defines the size of the final graph to be plotted for
%iterations
FixedStep = zeros(101,101);

for i = 1:100
    for j = 1:100
        X(:,:,1) = [x(i);y(j)];
        g(:,:,1) = [G{1}(x(i),y(j)),G{2}(x(i),y(j))];
        gT(:,:,1) = transpose(g(:,:,1));
        %the Q matrix for this function is [1,35;35,1);
        %the step size for alhpa that allows for convergence is
        %0<alpha<2/lambda(max)
        %lambda = 
        alpha = 1/2;
        for k = 2:100
            X(:,:,k) = X(:,:,k-1)+alpha*(-gT(:,:,k-1));
            g(:,:,k) = [G{1}(X(1,1,k),X(2,1,k)),G{2}(X(1,1,k),X(2,1,k))];
            gT(:,:,k) = transpose(g(:,:,k));
            if X(:,:,k)-X(:,:,k-1)<tolerance
                FixedStep(i,j) = k-1;
                break
            end
        end
    end
end

            




