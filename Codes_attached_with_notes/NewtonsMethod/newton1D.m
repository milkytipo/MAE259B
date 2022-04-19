% Filename: newton1D.m
% Demo of Newton-Raphson's method in 1D
%
% f(x) = x^2 - 3*x + 2 = 0
% f'(x) = 2*x - 3

function x = newton1D(  )

% Guess solution
x = 5;

% Tolerance
eps = 1e-6;

% Error
err = eps*10; % initialize to a value larger than eps

while err > eps
    deltaX = f(x) / df(x);
    x = x - deltaX;
    err = abs( f(x) );
end

fprintf('Solution is %f\n', x);

end

function s = f(x)
s = x^2 - 3*x + 2;
end

function s = df(x)
s = 2*x - 3;
end

