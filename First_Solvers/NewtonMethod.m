function root = NewtonsMethod(f, x0, tol, maxIter)
    for i = 1:maxIter
        fx = f(x0);
        dfx = (f(x0 + tol) - fx) / tol; % Numerical derivative
        x1 = x0 - fx / dfx;
        if abs(x1 - x0) < tol
            root = x1;
            return;
        end
        x0 = x1;
    end
    root = x0;
end