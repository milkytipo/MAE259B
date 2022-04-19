%% Global variables

%% Physical parameters
% Number of vertices
N = 3;

% Time step size
dt = 0.0000z1; % second

% Rod length
RodLength = 0.1; % meter

% Discrete length
deltaL = RodLength / (N-1);

% Radius of spheres
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;

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

% Total time
totalTime = 10; % seconds

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
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;

% Viscous damping matrix
C = zeros(6,6);
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Gravity
W = zeros(2*N,1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

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
all_mid_y = zeros( Nsteps, 1); % y-position of R2
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2

all_mid_y(1) = q(4);
all_mid_v(1) = u(4);

% Tolerance
tol = EI / RodLength^2 * 1e-3;

% Time marching scheme
for c=2:Nsteps
    
    fprintf('Time = %f\n', (c-1) * dt );
    
    q = q0; % Guess

    f = zeros(2*N,1);

    cofficient_q = M + C*dt;
    
    % Linear spring 1 between nodes 1 and 2
    xk = q(1);
    yk = q(2);
    xkp1 = q(3);
    ykp1 = q(4);
    f_Es_1 = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
    f(1:4) = f(1:4) + f_Es_1;

    % Linear spring 2 between nodes 2 and 3
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    f_Es_2 = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
    f(3:6) = f(3:6) + f_Es_2;
    
    % Bending spring between nodes 1, 2, and 3
    xkm1 = q(1);
    ykm1 = q(2);
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    curvature0 = 0;
    f_Eb = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
        curvature0, deltaL, EI);
    f(1:6) = f(1:6) + f_Eb;
    
    q = cofficient_q\(cofficient_q*q0 + M*u*dt - f*dt*dt + W*dt*dt);

    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
    
    figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis equal
    drawnow
    
    % Store
    all_mid_y(c) = q(4);
    all_mid_v(c) = u(4);
end

figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
