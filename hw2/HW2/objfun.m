function [x, d1, d2, refTwist] = objfun(x0, u, d1_old, d2_old, refTwist_old)

%% Global variables
global m dt unconsInd tol ScaleSolver maximum_iter
global Fg

iter = 0;
err = tol * ScaleSolver * 10;

x = x0; % Guess

while err > tol * ScaleSolver    
    iter = iter + 1;
    
    [d1, d2] = computeTimeParallel(d1_old, x0, x);
    tangent = computeTangent(x);
    refTwist = getRefTwist(d1, tangent, refTwist_old);
    theta = x(4:4:end);
    [m1,m2] = computeMaterialDirectors(d1, d2, theta);
    
    % Get forces
    [Fb, Jb] = getFb(x, m1, m2); % Bending
    [Ft, Jt] = getFt(x, refTwist); % Twisting
    [Fs, Js] = getFs(x);
    
    F = m .* (x-x0)/dt^2 - m .* u/dt - (Fb+Ft+Fs+Fg);    
    J = diag(m)/dt^2 - (Jb+Jt+Js);
    
    F_free = F(unconsInd);
    J_free = J(unconsInd, unconsInd);
    dx = J_free \ F_free;
    
    x(unconsInd) = x(unconsInd) - dx;
    
    err = sum(abs(F_free));
    fprintf('err = %f, iter= %d\n', err, iter);
    
    if iter > maximum_iter
        fprintf('Error\n');
        return
    end
end

end
