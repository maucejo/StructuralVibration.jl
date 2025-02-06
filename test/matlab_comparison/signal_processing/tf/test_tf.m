fs = 1024;
t = 0:1/fs:5-1/fs;
signal1 = randn(length(t), 1);
signal2 = cos(2*pi*10*t);
[H, f] = tfestimate(signal1, signal2, 512, 0, [], fs);

plot(f, abs(H))

% save test_tf.mat H f signal1 signal2