%this funciton will calcuate the hessian of the function f at the point (x0,y0)
function [H] = Hessian(f,x0,y0)
    syms x y 
    g = diff(f,x);
    gg = diff(g,x);
    h = diff(f,y);
    hh = diff(h,y);
    gh = diff(g,y);
    hg = diff(h,x);
    XX = gg(x0,y0);
    XY = gh(x0,y0);
    YX = hg(x0,y0);
    YY = hh(x0,y0);
    H = [XX,XY;YX,YY];
end