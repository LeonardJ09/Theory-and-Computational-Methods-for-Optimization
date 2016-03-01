%{
Jacob Leonard
MATH 467 - Fall 2015
jaleonar@usc.edu
Revision History
Date                    Changes                     Programmer
------------------------------------------------------------------
11/6/2015               Original                  Jacob Leonard
11/7/2015          Developed Derivatives          Jacob Leonard
11/10/2015        Developed Algorithm Body        Jacob Leonard
11/12/2015        Figured Out Alpha Threshold     Jacob Leonard
11/13/2015          Troubleshooting               Jacob Leonard
11/14/2015          Developed Z Function          Jacob Leonard

%}

%this script is for the fixed step size gradient method

%determine x(0) and y(0) for the start of the methods
for i = 1:101
    x(i) = (-2)+((4*(i-1))/100);
    y(i) = (-2)+((4*(i-1))/100);
end

%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(x,y) ((x^4+y^4-6*x^2*y^2-1)^2+(4*x^3*y-4*x*y^3)^2);
G = {@(x,y) (8*x*(x^6+3*x^4*y^2+x^2*(3*y^4-1)+y^2*(y^4+3))),@(x,y) (8*y*(x^6+3*x^4*y^2+3*x^2*(y^4+1)+y^2*(y^4-1)))};
%Gradient = [g{1}(x,y),g{2}(x,y)];
%when the
H = {@(x,y) 8*(7*x^6+15*x^4*y^2+x^2*(9*y^4-3)+y^2*(y^4+3)),@(x,y) 48*x*y*(x^4+2*x^2*y^2+y^4+1);@(x,y) 48*x*y*(x^4+2*x^2*y^2+y^4+1),@(x,y) 8*(x^6+9*x^4*y^2+3*x^2*(5*y^4+1)+y^2*(7*y^4-3))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%desired level of accuracy
tolerance = 10^(-7);

%this matrix defines the size of the final graph to be plotted for
%iterations
FixedStep = zeros(101,101);

for i = 1:101
    for j = 1:101
        X(:,:,1) = [x(i);y(j)];
        g(:,:,1) = [G{1}(x(i),y(j)),G{2}(x(i),y(j))];
        gT(:,:,1) = transpose(g(:,:,1));
        %the Q matrix for this function is [1,35;35,1);
        alpha = .1;
        for k = 2:100
            X(:,:,k) = X(:,:,k-1)-alpha*(gT(:,:,k-1));
            g(:,:,k) = [G{1}(X(1,1,k),X(2,1,k)),G{2}(X(1,1,k),X(2,1,k))];
            gT(:,:,k) = transpose(g(:,:,k));
            if f(X(1,1,k),X(2,1,k))-f(X(1,1,k),X(2,1,k))<tolerance
                FixedStep(i,j) = k;
                break
            end
        end
    end
end

            




