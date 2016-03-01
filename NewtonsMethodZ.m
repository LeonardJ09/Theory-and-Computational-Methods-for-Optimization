%{
Jacob Leonard
MATH 467 - Fall 2015
jaleonar@usc.edu
Revision History
Date             Changes                  Programmer
-----------------------------------------------------
11/15/2015        Original                Jacob Leonard
       ~Seeking Counsel from Professor Wang~
12/4/2015      Reverted Function          Jacob Leonard
12/4-12/2015    Troubleshooting           Jacob Leonard
12/13/2015    Fixing Backtracking         Jacob Leonard
12/15/2015 Analyzing Backtracking Values  Jacob Leonard      
12/16/2015        Completed               Jacob Leonard
%}

%this script is for newtons method with back-tracking

%define the values of x and y from -2 to 2, increasing by 1/25, for 101
%values
for j = 1:101
    x(j) = (-2+(4*(j-1))/100);
    y(j) = (-2+(4*(j-1))/100);
end

%define an anonymous function handle for the equations that compose the gradient and the hessian
f = @(x,y) ((x^4+y^4-6*x^2*y^2-1)^2+(4*x^3*y-4*x*y^3)^2);
G = {@(x,y) (8*x*(x^6+3*x^4*y^2+x^2*(3*y^4-1)+y^2*(y^4+3))),@(x,y) (8*y*(x^6+3*x^4*y^2+3*x^2*(y^4+1)+y^2*(y^4-1)))};
%Gradient = [g{1}(x,y),g{2}(x,y)];
H = {@(x,y) (8*(7*x^6+15*x^4*y^2+x^2*(9*y^4-3)+y^2*(y^4+3))),@(x,y) (48*x*y*(x^4+2*x^2*y^2+y^4+1));@(x,y) (48*x*y*(x^4+2*x^2*y^2+y^4+1)),@(x,y) (8*(x^6+9*x^4*y^2+3*x^2*(5*y^4+1)+y^2*(7*y^4-3)))};
%Hessian = [H{1}(x,y),H{2}(x,y);H{3}(x,y),H{4}(x,y)];

%desired level of accuracy
tolerance = 10^(-7);

%this matrix shows the number of iterations for matlab to think the value
%is zero, or within the desired tolerance
NewtonsSteps = zeros(101,101);
%this is the value of the function at the point the algorithm terminated
NewtonsValues = zeros(101,101);


%given each initial value, newtons method will iterate without the need to
%evaluate any points

for i = 1:101
    for j = 1:101
        Z(:,:,1) = [x(i);y(j)];
        %set the values in Z equal to 0 to track the progress
        Z(1,1,2:5000)=0;
        Z(2,1,2:5000)=0;
        g(:,:,1) = [G{1}(x(i),y(j)),G{2}(x(i),y(j))];
        %if the gradient is equal to zero, then we have reached the
        %optimal solution according to the algorithm
        if (g(1,1,1) == 0) && (g(1,2,1) == 0)
            NewtonsSteps(i,j) = 0;
            NewtonsValues(i,j) = f(Z(1,1,1),Z(2,1,1));
            continue
        end
        gT(:,:,1) = transpose(g(:,:,1));
        h(:,:,1) = [H{1}(x(i),y(j)),H{2}(x(i),y(j));H{3}(x(i),y(j)),H{4}(x(i),y(j))];
        I(:,:,1) = pinv(h(:,:,1));
        d(:,:,1)=(I(:,:,1)*gT(:,:,1));
        %backtracking values
        B = .5;
        A = .1;
        t1(1)=1;
        for k = 2:5000
            %add backtracking to the line search
            s = 0;
            %return the backtracking values to 0
            t1(k) = 1;
            while s==0
                %in order for the backtracking search to work, the
                %following rule must be satisfied: m>=c
                c = f(Z(1,1,k-1)-t1(k)*d(1,1,k-1),Z(2,1,k-1)-t1(k)*d(2,1,k-1));
                m = f(Z(1,1,k-1),(Z(2,1,k-1)))-A*t1(k)*g(:,:,k-1)*d(:,:,k-1);
                v = c-m;
                if v <= tolerance
                    t2 = t1(k);
                    break
                end
                t1(k) = B*t1(k);
            end
            %find the new Z value with the updated t(k)
            Z(:,:,k) = Z(:,:,k-1)-t2*d(:,:,k-1);
            %check to see if the Z values have NaN or Inf
            %values
            if (isnan(Z(1,1,k)) == 1) || (isnan(Z(2,1,k)) == 1)
                NewtonsSteps(i,j) = 5000;
                NewtonsValues(i,j) = 1;                                                                            ;
                break
            end
            if (isinf(Z(1,1,k)) == 1) || (isinf(Z(2,1,k)) == 1)
                NewtonsSteps(i,j) = 5000;
                NewtonsValues(i,j) = 1;                                                                            ;
                break
            end
            %if the function value dips below the tolerance, then it is
            %considered to have converged to the optimal value
            if f(Z(1,1,k),Z(2,1,k))<tolerance;
                NewtonsSteps(i,j) = k-1;
                NewtonsValues(i,j) = 0;
                break
            end
            %calculate the new gradient at the updated point
            g(:,:,k) = [G{1}(Z(1,1,k),Z(2,1,k)),G{2}(Z(1,1,k),Z(2,1,k))];
            %if the gradient is equal to zero, then we have reached the
            %optimal solution according to the algorithm
            if (g(1,1,k)) == 0 && (g(1,2,k) == 0)
                NewtonsSteps(i,j) = k-1;
                NewtonsValues(i,j) = f(Z(1,1,k),Z(2,1,k));                                                                            ;
                break
            end
            %check to see if the gradient has NaN or Inf
            %values
            if (isnan(g(1,1,k)) == 1) || (isnan(g(1,2,k)) == 1)
                NewtonsSteps(i,j) = 5000;
                NewtonsValues(i,j) = 1;                                                                            ;
                break
            end
            if (isinf(g(1,1,k)) == 1) || (isinf(g(1,2,k)) == 1)
                NewtonsSteps(i,j) = 5000;
                NewtonsValues(i,j) = 1;                                                                            ;
                break
            end
            h(:,:,k) = [H{1}(Z(1,1,k),Z(2,1,k)),H{2}(Z(1,1,k),Z(2,1,k));H{3}(Z(1,1,k),Z(2,1,k)),H{4}(Z(1,1,k),Z(2,1,k))];
            I(:,:,k) = pinv(h(:,:,k));
            gT(:,:,k) = transpose(g(:,:,k));
            d(:,:,k)=(I(:,:,k)*gT(:,:,k));
            if k == 5000
                NewtonsSteps(i,j) = 5000;
                NewtonsValues(i,j) = 1;
                break
            end
        end
    end
end

xAxis = linspace(-2,2,101);
yAxis = linspace(-2,2,101);

%this plot will look at the number of steps it took for the algorithm to finish, the real values of the function over the interval
%for which the algorithm finished, and the  imaginary values for which the
%algorithm finished

subplot(2,2,1:2)
%this plot will show the number of iterations it took
contourf(xAxis,yAxis,NewtonsSteps);
xlabel('x');
ylabel('y');
title('Newtons Method with Backtracking # of Steps, B=.5, A =.1');
colorbar;
subplot(2,2,3:4)
NewtonsValuesReal = real(NewtonsValues);
contourf(xAxis,yAxis,NewtonsValuesReal);
xlabel('x');
ylabel('y');
title('Binary Convergence Plot x=[-2:2], y=[-2:2]');
colorbar;

