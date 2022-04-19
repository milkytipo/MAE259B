%% Global variables

%% Physical parameters
function v_terminal = Eleven_Nodes_Implicit(N, dt)

plot_subfig = 1;

if plot_subfig == 1
    N = 21;
    dt = 0.1; % second
end

% Total time
totalTime = 50; % seconds

% Rod length
RodLength = 0.1; % meter

% Discrete length
deltaL = RodLength / (N-1);

% Radius of spheres
R = ones(N,1);
R = R*deltaL/10;
R(round(N/2)) = 0.025;

% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

% Rod radius
r0 = 0.001;

% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Viscosity
visc = 1000; % Pa-s

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;

% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
%     nodes(c,2) = 0;
end

% Mass matrix
M = zeros(2*N,2*N);

for i = 1:2:2*N
    j = round(i/2);
    M(i,i) = 4/3*pi*R(j)^3*rho_metal;
    M(i+1,i+1) = 4/3*pi*R(j)^3*rho_metal;
end

% Viscous damping matrix
C = zeros(2*N,2*N);

for i = 1:2:2*N
    j = round(i/2);
    C(i,i) = 6*pi*visc*R(j);
    C(i+1,i+1) = 6*pi*visc*R(j);
end

% Gravity
W = zeros(2*N,1);

for i = 1:N
    j=2*i;
    W(j) = -4/3*pi*R(i)^3*rho*g;
end

% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R_mid
all_mid_v = zeros( Nsteps, 1); % y-velocity of R_mid

all_mid_y(1) = q(2*round(N/2));
all_mid_v(1) = u(2*round(N/2));

% Tolerance
tol = EI / RodLength^2 * 1e-3;

% Time marching scheme
for c=2:Nsteps
    
%     fprintf('Time = %f\n', (c-1) * dt );
    
    q = q0; % Guess
    
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        % Inertia
        f = M / dt * ( (q-q0) / dt - u );
        J = M / dt^2;
        
        %
        % Elastic forces
        %
        % Linear spring force
        for i = 1:2:2*(N-1) 
            xk = q(i);
            yk = q(i+1);
            xkp1 = q(i+2);
            ykp1 = q(i+3);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            f(i:i+3) = f(i:i+3) + dF;
            J(i:i+3,i:i+3) = J(i:i+3,i:i+3) + dJ;         
        end
        
        % Bending spring between nodes 1, 2, and 3
        for i = 1:2:2*(N-2)
            xkm1 = q(i);
            ykm1 = q(i+1);
            xk = q(i+2);
            yk = q(i+3);
            xkp1 = q(i+4);
            ykp1 = q(i+5);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                curvature0, deltaL, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                curvature0, deltaL, EI);
            f(i:i+5) = f(i:i+5) + dF;
            J(i:i+5,i:i+5) = J(i:i+5,i:i+5) + dJ;
        end
        
        % Viscous force
        f = f + C * ( q - q0 ) / dt;
        J = J + C / dt;
        
        % Weight
        f = f - W;
        
        % Update
        q = q - J \ f;
        
        err = sum( abs(f) );
    end

    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
    
    if plot_subfig == 1
        figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        drawnow
    end
    % Store
    all_mid_y(c) = q(2*round(N/2));
    all_mid_v(c) = u(2*round(N/2));
end
v_terminal = all_mid_v(end)

if plot_subfig == 1
    figure(2);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_v, 'k-');
    xlabel('Time, t [sec]');
    ylabel('Velocity of mid-node, v [meter/sec]');

    
    figure(3);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_y, 'k-');
    xlabel('Time, t [sec]');
    ylabel('Position of mid-node, y [meter]');
end
end
