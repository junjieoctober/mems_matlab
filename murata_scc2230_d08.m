%analyze murata acc
%the file is around 50minutes with the murata at the corner of zhai's desk
function murata_scc2230_d08()

clc;

%close all; 


%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_50min_dec4_2015.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_2min_dec4_2015_1stcw_2ndccw.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_free_drop.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_overnight_10hrs_dec4_2015.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_2hrs_dec5_2015.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_3hrs_dec6_2015.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_walking.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec7_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec8_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec9_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec10_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec11_2days.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec15_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec16_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec17_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec18_3days.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
%[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec22_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);
[time,celsius,x_gyro,y_gyro,z_gyro, x_acc,y_acc,z_acc] = textread('data/capture_dec23_overnight.txt', '$DRAW,%f,%f,%f,%f,%f,%f,%f,%f',-1);

%On Time (s),Algorithm,Coordinates,Roll (deg),Pitch (deg),Compass (deg),Accelerometer (g) X,Accelerometer (g) Y,Accelerometer (g) Z,Magnetometer (uT) X,Magnetometer (uT) Y,Magnetometer (uT) Z,Gyroscope (deg/s) X,Gyroscope (deg/s) Y,Gyroscope (deg/s) Z,q0,q1,q2,q3,Ang Vel (deg/s) X,Ang Vel (deg/s) Y,Ang Vel (deg/s) Z,Acceleration (g) X,Acceleration (g) Y,Acceleration (g) Z,Hard iron (uT) X,Hard iron (uT) Y,Hard iron (uT) Z,Soft iron XX,Soft iron YY,Soft iron ZZ,Soft iron XY,Soft iron XZ,Soft iron YZ,Measurements,Fit Error %,Geomagnetic field (uT),Inclination Delta (deg),q_zge X,q_zge Y,q_zge Z,q_zme X,q_zme Y,q_zme Z,q_ge+ X,q_ge+ Y,q_ge+ Z,q_me+ X,q_me+ Y,q_me+ Z,Offset b+ (deg/s) X,Offset b+ (deg/s) Y,Offset b+ (deg/s) Z,Delta+ (deg),Altitude (m),Temp (C)
%[time,roll,pitch,compass,x_acc,y_acc,z_acc,MagnetometerX,MagnetometerY,MagnetometerZ,x_gyro,y_gyro,z_gyro,q0,q1,q2,q3,AngVelX,AngVelY,AngVelZ,AccelerationX,AccelerationY,AccelerationZ,HardironX,HardironY,HardironZ,SoftironXX,SoftironYY,SoftironZZ,SoftironXY,SoftironXZ,SoftironYZ,Measurements,FitError,Geomagneticfield,InclinationDelta,q_zgeX,q_zgeY,q_zgeZ,q_zmeX,q_zmeY,q_zmeZ,q_geX,q_geY,q_geZ,q_meX,q_meY,q_meZ,Offset_bX,Offset_bY,Offset_bZ,Delta,Altitude,celsius] = textread('data/sensorfusion_dec23_overnight_s01.csv', '%f,Acc+Mag+Gyro,Aerospace,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f','headerlines', 1);
%[time,roll,pitch,compass,x_acc,y_acc,z_acc,MagnetometerX,MagnetometerY,MagnetometerZ,x_gyro,y_gyro,z_gyro,q0,q1,q2,q3,AngVelX,AngVelY,AngVelZ,AccelerationX,AccelerationY,AccelerationZ,HardironX,HardironY,HardironZ,SoftironXX,SoftironYY,SoftironZZ,SoftironXY,SoftironXZ,SoftironYZ,Measurements,FitError,Geomagneticfield,InclinationDelta,q_zgeX,q_zgeY,q_zgeZ,q_zmeX,q_zmeY,q_zmeZ,q_geX,q_geY,q_geZ,q_meX,q_meY,q_meZ,Offset_bX,Offset_bY,Offset_bZ,Delta,Altitude,celsius] = textread('data/sensorfusion_dec23_overnight.csv', '%f,Acc+Mag+Gyro,Aerospace,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f','headerlines', 1);
%[time,roll,pitch,compass,x_acc,y_acc,z_acc,MagnetometerX,MagnetometerY,MagnetometerZ,x_gyro,y_gyro,z_gyro,q0,q1,q2,q3,AngVelX,AngVelY,AngVelZ,AccelerationX,AccelerationY,AccelerationZ,HardironX,HardironY,HardironZ,SoftironXX,SoftironYY,SoftironZZ,SoftironXY,SoftironXZ,SoftironYZ,Measurements,FitError,Geomagneticfield,InclinationDelta,q_zgeX,q_zgeY,q_zgeZ,q_zmeX,q_zmeY,q_zmeZ,q_geX,q_geY,q_geZ,q_meX,q_meY,q_meZ,Offset_bX,Offset_bY,Offset_bZ,Delta,Altitude,celsius] = textread('data/sensorfusion_dec28_overnight.csv', '%f,Acc+Mag+Gyro,Aerospace,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f','headerlines', 1);
%[time,roll,pitch,compass,x_acc,y_acc,z_acc,MagnetometerX,MagnetometerY,MagnetometerZ,x_gyro,y_gyro,z_gyro,q0,q1,q2,q3,AngVelX,AngVelY,AngVelZ,AccelerationX,AccelerationY,AccelerationZ,HardironX,HardironY,HardironZ,SoftironXX,SoftironYY,SoftironZZ,SoftironXY,SoftironXZ,SoftironYZ,Measurements,FitError,Geomagneticfield,InclinationDelta,q_zgeX,q_zgeY,q_zgeZ,q_zmeX,q_zmeY,q_zmeZ,q_geX,q_geY,q_geZ,q_meX,q_meY,q_meZ,Offset_bX,Offset_bY,Offset_bZ,Delta,Altitude,celsius] = textread('data/sensorfusion_dec29_overnight.csv', '%f,Acc+Mag+Gyro,Aerospace,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f','headerlines', 1);



figure;
plot(celsius);
title('temperature'); 


figure;
plot(x_gyro);
hold on;
plot(y_gyro);
plot(z_gyro); 


x_g_lpf = zeros(length(x_gyro),1); 
y_g_lpf = zeros(length(y_gyro),1); 
z_g_lpf = zeros(length(z_gyro),1); 

x_g_lpf(1) = x_gyro(1); 
y_g_lpf(1) = y_gyro(1); 
z_g_lpf(1) = z_gyro(1); 

lpf = 0.001; 
for i = 2: length(x_gyro)
    x_g_lpf(i) = x_g_lpf(i-1) + lpf*(x_gyro(i)-x_g_lpf(i-1));
    y_g_lpf(i) = y_g_lpf(i-1) + lpf*(y_gyro(i)-y_g_lpf(i-1));
    z_g_lpf(i) = z_g_lpf(i-1) + lpf*(z_gyro(i)-z_g_lpf(i-1));
end
plot(x_g_lpf);
hold on;
plot(y_g_lpf);
plot(z_g_lpf); 
%ylim([-3e-3 3e-3]); 

title('gyro raw');
legend('x','y','z');
str=sprintf('mx\\_gyro = %f, stdx\\_gyro = %f,my\\_gyro = %f,\n stdy\\_gyro = %f,mz\\_gyro = %f,stdz\\_gyro = %f',mean(x_gyro),std(x_gyro),mean(y_gyro),std(y_gyro),mean(z_gyro),std(z_gyro)); 
display(str); 
gtext(str); 


% h = 0.1; 
% grad_x_gyro = gradient(x_gyro,h); 
% m_grad_x_gyro = mean(grad_x_gyro);
% grad_y_gyro = gradient(y_gyro,h); 
% m_grad_y_gyro = mean(grad_y_gyro);
% grad_z_gyro = gradient(z_gyro,h); 
% m_grad_z_gyro = mean(grad_z_gyro);
% str=sprintf('m_grad_x_gyro = %f,m_grad_y_gyro = %f,m_grad_z_gyro = %f',m_grad_x_gyro,m_grad_y_gyro,m_grad_z_gyro); 
% display(str); 
% gtext(str); 

% 
% idx1 = find(time==19481.05); 
% idx2 = find(time==19492.05);
% v = (idx1:idx2); 
% s = trapz(time(v),z_gyro(v)-(-0.022940)); 
% str=sprintf('integral angle = %f',s); 
% display(str); 
% 
% idx1 = find(time==19511.05); 
% idx2 = find(time==19529.05);
% v = (idx1:idx2); 
% s = trapz(time(v),z_gyro(v)-(-0.022940)); 
% str=sprintf('integral angle = %f',s); 
% display(str); 


x_a_lpf = zeros(length(x_acc),1); 
y_a_lpf = zeros(length(y_acc),1); 
z_a_lpf = zeros(length(z_acc),1); 
x_a_lpf(1) = x_acc(1); 
y_a_lpf(1) = y_acc(1); 
z_a_lpf(1) = z_acc(1); 

for i = 2: length(x_acc)
    x_a_lpf(i) = x_a_lpf(i-1) + lpf*(x_acc(i)-x_a_lpf(i-1));
    y_a_lpf(i) = y_a_lpf(i-1) + lpf*(y_acc(i)-y_a_lpf(i-1));
    z_a_lpf(i) = z_a_lpf(i-1) + lpf*(z_acc(i)-z_a_lpf(i-1));
end

figure;
subplot(3,1,1);
plot(x_acc);
hold on; 
plot(x_a_lpf);



subplot(3,1,2); 
plot(y_acc);
hold on; 
plot(y_a_lpf);


subplot(3,1,3); 
plot(z_acc);
hold on; 
plot(z_a_lpf);


title('acc raw'); 
legend('x','y','z');
str=sprintf('mx\\_acc = %f, stdx\\_acc = %f,my\\_acc = %f,\n stdy\\_acc = %f,mz\\_acc = %f,stdz\\_acc = %f',mean(x_acc),std(x_acc),mean(y_acc),std(y_acc),mean(z_acc),std(z_acc)); 
display(str); 
gtext(str); 

dt=time(end)-time(1); 
str=sprintf(' duration = %f seconds',dt);
display(str); 




%next calculate Allan-Variance
T0 = time(2)- time(1); 
fs = 1/T0; 
pts = 100; 
[T,sigma] = allan(z_gyro,fs,pts)
figure; 
loglog(T,sigma,'*-'); 
grid on; 


figure; 
timezero = (1:ceil(60/T0));
z_gyrorw = z_gyro - mean(z_gyro(timezero)); 
y_gyrorw = y_gyro - mean(y_gyro(timezero)); 
x_gyrorw = x_gyro - mean(x_gyro(timezero)); 

thetaz= cumsum(z_gyrorw*T0);
thetay= cumsum(y_gyrorw*T0);
thetax= cumsum(x_gyrorw*T0);
plot(thetaz); 
grid on;
hold on;
plot(thetay);
plot(thetax);
title('z gyro angle ');
legend('z gyro','y gyro','x gyro'); 