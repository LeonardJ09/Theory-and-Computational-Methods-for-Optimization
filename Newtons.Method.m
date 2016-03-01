%{ 
Jacob Leonard
MATH 467 - Fall 2015
jaleonar@usc.edu
Revision History
Date             Changes           Programmer 
-------------------------------------------------
11/6/2015        Original          Jacob Leonard

%}

%this script is for newtons method with back-tracking

%determine x(0) and y(0) for the start of the methods
for j = 1:100
    x(j) = (-2)+((4*j)/100);
    y(j) = (-2)+((4*j)/100);
end


%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(x,y) (abs(x^4 + y^4 -6*x^2*y^2-1)^2+1i*(4*x^3*y-4*x*y^3)-1)^2;
G = {@(x,y) (8*((x+1i*y)^4-1)*(x+1i*y)^2);@(x,y) (8*1i*((x+1i*y)^4-1)*(x+1i*y)^2)};
%Gradient = [g{1}(x,y),g{2}(x,y)];
%when the 
H = {@(x,y) (24*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*1i*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*(x+1i*y)^2*((7/3)*(x+1i*y)^4-1));@(x,y) (24*(x+1i*y)^2*(-(7/3)*(x+1i*y)^4+1))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%this matrix defines the size of the final graph to be plotted for
%iterations
Newtons = zeros(100,100);

%desired level of accuracy
tolerance = 10^(-7);

%lowest value we wish to divide by
epsilon = 10^(-14);

%given each initial value, newtons method will iterate without the need to
%evaluate any points
for i = 1:100
    for j = 1:100
            X(:,:,1) = [x(i);y(j)];
            g(:,:,1) = [G{1}(x(i),y(j)),G{2}(x(i),y(j))];
            h(:,:,1) = [H{1}(x(i),y(j)),H{2}(x(i),y(j));H{3}(x(i),y(j)),H{4}(x(i),y(j))]; 
            for k = 2:100
                X(:,:,k) = X(:,:,k-1)-(((g(:,:,k-1)/(h(:,:,k-1))))));
                g(:,:,k) = [G{1}(X(1,1,k),X(2,1,k))),G{2}(X(1,1,k),X(2,1,k))];
                h(:,:,k) = [H{1}((X(1,1,k),X(2,1,k)),H{2}((X(1,1,k),X(2,1,k));H{3}((X(1,1,k),X(2,1,k)),H{4}((X(1,1,k),X(2,1,k))];
                if X(:,:,k)-X(:,:,k-1) < tolerance
                    Newtons(i,j) = k;
                    break
                end
                if h(:,:,k) < epsilon
                    Newtons(i,j) = k;
                    break
                end
            end
    end
end
                