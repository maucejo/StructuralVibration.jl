fs = 2048;
fn = fs/2;

% Filter design
order = 200;
fb = 0.05*fn;
b = firls(order, [0 1/2-0.05 1/2+0.05 1] , [1, 1, 0, 0]);

[h, w] = freqz(b, 1);
freq = w/pi*fs;

% Visualization
figure
subplot(2, 1, 1)
plot(freq, abs(h))
subplot(2, 1, 2)
plot(freq, unwrap(angle(h)))

hf = h(freq <= fn);
freq = freq(freq <= fn);

% save test_anti_aliasing.mat hf freq