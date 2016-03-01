%this funciton will calcuate the hessian of the function f at the point (x0,y0)
function [D] = FirstOrder(f,x0,y0)
    syms x y 
    g = diff(f,x);
    h = diff(f,y);
    X = g(x0,y0);
    Y = h(x0,y0);
    D = [X,Y];
end