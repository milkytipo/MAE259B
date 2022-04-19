
seq_p = (500:1000:20000);

y_simu = zeros(size(seq_p));
y_pred = zeros(size(seq_p));
index_P = 1;

for P = seq_p
    [y_simu(index_P), y_pred(index_P)] = euler_bernoulli_beam(P);
    index_P = index_P+1;
end

figure(1);
plot(seq_p, y_simu, 'k-');
hold on;
plot(seq_p, y_pred, 'g-');
xlabel('P[N]'); 

ylabel('Comparison between y_red and y_simu, v [meter/sec]');