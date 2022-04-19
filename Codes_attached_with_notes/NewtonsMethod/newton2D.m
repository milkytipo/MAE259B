% Filename: newton2D.m
% Demo of Newton-Raphson's method in 2D
%
% f(x) = [   x1^3+ x2   - 1;
%           -x1  + x2^3 + 1  ]
%
% J(x) = [  3x1^2,  1;
%           -1,     -3x2^2  ]
%

function x = newton2D(  )

% Guess solution
x = [0.5; 0.5];

% Tolerance
eps = 1e-6;

% Error
err = eps*10; % initialize to a value larger than eps

while err > eps
    deltaX = J(x) \ f(x);
    x = x - deltaX;
    err = abs( sum( f(x) ) );
end

fprintf('Solution is x1=%f, x2=%f\n', x(1), x(2));

end

function s = f(x)
s = [x(1)^3 + x(2) - 1; ...
    -x(1) + x(2)^3 + 1];
end

function s = J(x)
s = [   3*x(1)^2,       1;      ...
        -1,             3*x(2)^2];
end

