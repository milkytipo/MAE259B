%% Global variables

%% Physical parameters
% default: Number of vertices and Time step size
% N = 50;
% 
% 
% % terminal speed vs. dt
% M = 5;
% v_terminal = zeros(M,1);
% index_M = 1;

% seq_dt = (0.0001:0.005:0.1);
% for dt = seq_dt
%     v_terminal(index_M) = Eleven_Nodes_Implicit(N, dt);
%     index_M = index_M+1;
% end
% 
% figure(1);
% plot(seq_dt, v_terminal, 'k-');
% xlabel('dt, t [sec]'); 
% ylabel('Velocity of mid-node, v [meter/sec]');

% terminal speed vs. node numbers
v_terminal = zeros(N,1);
dt = 0.01; % second

for N = 3:21
    v_terminal(N) = Eleven_Nodes_Implicit(N, dt);
end
figure(2);
timeArray = 3:21;
plot(timeArray, v_terminal(3:21), 'k-');
xlabel('Node, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');

