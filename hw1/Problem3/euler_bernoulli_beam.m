%% Global variables

function [y_simu, y_pred] = euler_bernoulli_beam(P)
plot_subfig = 0;

if plot_subfig == 1
    P = 20000;
end

% Total time
totalTime = 0.5; % seconds

N = 50;
dt = 0.0005; % second
% Rod length
RodLength = 1; % meter
R = 0.013;
r = 0.011;

P_position = 0.75;
E = 7e10;
I = pi/4.0*(R.^4-r.^4);
rou = 2700;
m = pi*(R*R-r*r)*RodLength*rou/(N-1);


% Discrete length
deltaL = RodLength / (N-1);
index_node_P = round(P_position/deltaL);

% Young's modulus
Y = E; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Viscosity
visc = 1000; % Pa-s

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * I;
EA = Y * pi*(R*R-r*r);

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
    M(i,i) = m;
    M(i+1,i+1) = m;
end

% Viscous damping matrix
C = zeros(2*N,2*N);

% Gravity
W = zeros(2*N,1);

for i = 1:N
    j=2*i;
    W(j) = m*g;
end
W(2*index_node_P) = -P;

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
    
    fprintf('Time = %f\n', (c-1) * dt );
    
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
            if i==1
                f(3:i+3) = f(3:i+3) + dF(3:4);
                J(3:i+3,3:i+3) = J(3:i+3,3:i+3) + dJ(3:4,3:4);  
            elseif i == 2*(N-1)-1
                f(i:i+2) = f(i:i+2) + dF(1:3);
                J(i:i+2,i:i+2) = J(i:i+2,i:i+2) + dJ(1:3,1:3);  
            else
                f(i:i+3) = f(i:i+3) + dF;
                J(i:i+3,i:i+3) = J(i:i+3,i:i+3) + dJ;  
            end
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

            if i==1
                f(3:6)= f(3:6) + dF(3:6);
                J(3:6,3:6) = J(3:6,3:6) + dJ(3:6,3:6);  
            elseif i == 2*(N-2)-1
                f(i:i+4) = f(i:i+4) + dF(1:5);
                J(i:i+4,i:i+4) = J(i:i+4,i:i+4) + dJ(1:5,1:5); 
            else
                f(i:i+5) = f(i:i+5) + dF;
                J(i:i+5,i:i+5) = J(i:i+5,i:i+5) + dJ;
            end           
        end
        
        % Weight
        f = f - W;
        
        % Update
%         q(3:2*N-1) = q(3:2*N-1) - J(3:2*N-1,3:2*N-1) \ f(3:2*N-1);
        q = q- J\ f;
        
        err = sum( abs(f) );
    end
    q(1) = 0;
    q(2) = 0;
    q(end) = 0;
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
v_terminal = all_mid_v(end);

if plot_subfig == 1
    figure(2);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_v, 'k-');
    xlabel('Time, t [sec]');
    ylabel('Velocity of mid-node, v [meter/sec]');
end

c_min = RodLength-P_position;
y_max = (P * c_min * (RodLength*RodLength-c_min*c_min).^1.5)/(9*sqrt(3)*E*I*RodLength);

y_simu = all_mid_y(end);
y_pred = y_max;

end
