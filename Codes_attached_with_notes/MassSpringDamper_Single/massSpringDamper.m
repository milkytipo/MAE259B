% Filename: massSpringDamper.m
% Simulation of a mass-spring-damper system
function massSpringDamper()
global m K b F0 omega dt

m = 1; % mass
K = 1; % spring constant
b = 10; % viscous damping coefficient
F0 = 0.1; % driving amplitude
omega = 1; % driving frequency
maxTime = 50; % total time of simulation
dt = 1e-1; % discrete time step
eps = 1e-6 * F0; % error tolerance

t = 0:dt:maxTime; % time
x = zeros(size(t)); % position
u = zeros(size(t)); % velocity

% Initial condition: position and velocity are zero
x(1) = 0;
u(1) = 0;

for k=1:length(t)-1 % march over time steps
    x_old = x(k);
    u_old = u(k);    
    t_new = t(k+1);
    
    x_new = x_old; % guess solution
    err = eps * 100; % initialize to a large value
    while err > eps
        deltaX = f(x_new, x_old, u_old, t_new) / ...
            df(x_new, x_old, u_old, t_new);
        x_new = x_new - deltaX;
        err = abs( f(x_new, x_old, u_old, t_new) );
    end
    
    u_new = (x_new - x_old) / dt;
    
    x(k+1) = x_new; % store solution
    u(k+1) = u_new;
end

% Plot it
plot(t, x, 'ro-'); xlabel('t [second]'); ylabel('x [meter]');
end

function s = f(x_new, x_old, u_old, t_new)
global m K b F0 omega dt
s = m/dt * ( (x_new-x_old)/dt - u_old ) + K*x_new + b*(x_new-x_old)/dt ...
    - F0*sin(omega*t_new);
end

function s = df(x_new, x_old, u_old, t_new)
global m K b dt
s = m/dt^2 + K + b/dt;
end
