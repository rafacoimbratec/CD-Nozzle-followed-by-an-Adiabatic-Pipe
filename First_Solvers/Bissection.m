function [x] = Bissection(f,a,b,n)
format long; %Use as many decimal places as possible
c = f(a); %Set c equal to function at a
d = f(b); %Set d equal to function at b
if c*d > 0.0 %If c times d is positive, bisection method will not work.
    fprintf('\nError\n');
    error('Function has the same sign at both endpoints.Bisection method will not work');
end
for k = 1:n %Perform n iterations of bisection method
    x = (a+b)/2; %new x is half the interval
    y = f(x); %Evaluate f(x)
    %disp([x y]); %display value of (x,y) coordinate
    if y == 0.0 %If y equals zero, have found root 
        e = 0; %Set initial error equal to zero
        break %Exit for loop
    end 
    if c*y<0 %Opposite sign
        b = x;%Set new right limit to x ((a+b)/2)
    else a = x;%Otherwise, set left limit to x (a=x);
    end
end
end