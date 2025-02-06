fs = 1024;
t = 0:1/fs:5-1/fs;
signal = cos(2*pi*100*t);

[pxx_mat, fmat] = pwelch(signal, 512, 256, 512, fs, 'power');

plot(fmat, 10*log10(pxx_mat))

save test_spectrum.mat pxx_mat fmat signal