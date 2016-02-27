function simulate_gyro_arw()

N = 1e5;

gyro = zeros(N,1);
gyro_lpf = zeros(N,1); 
np   = 1; 
bias = 1; 
Ts = 0.01; 
k = 0*Ts; 
time = Ts*(1:N); 
whitenoise = np*randn(N,1); 
rwnoise    = zeros(N,1); 
gyro(1) = bias;
gyro_lpf(1) = gyro(1); 
for i = 2 : N
    rwnoise(i) =  whitenoise(i); 
    gyro(i) = bias + k*i+ rwnoise(i);
    gyro_lpf(i) = gyro_lpf(i-1)+0.001*(gyro(i)-gyro_lpf(i-1)); 
end

%next calculate Allan-Variance
theta = cumsum(gyro); 
figure; 
plot(time,gyro); 
hold on;
plot(time,theta); 
plot(time,gyro_lpf); 
fs = 1/Ts; 
pts = 1000; 
[T,sigma] = allan(gyro,fs,pts)
figure; 
loglog(T,sigma,'*-'); 
grid on;