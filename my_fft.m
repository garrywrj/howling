x =audioread('howl.wav');

fs =16000;
N = length(x);
t = (0:N-1)/fs;
f = (0:N-1)*fs/N;

y = abs(fft(x,N));

f=f(1:N/2);
y=y(1:N/2);
plot(f,y);
