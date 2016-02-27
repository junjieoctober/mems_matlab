function analog_device_gyro_filter_simulate()
%http://www.analog.com/library/analogDialogue/archives/46-07/MEMS_stabilization.pdf

Fmax = 9840/2; % one-half of the sample rate
w = zeros(Fmax,1); 
for f = 1:Fmax
 w(f) = 2*pi*f;
end
p1 = 404; % pole location = 404Hz
p2 = 757; % pole location = 757Hz
NUM1 = 2*pi*p1;
DEN1 = [1 2*pi*p1];
NUM2 = 2*pi*p2;
DEN2 = [1 2*pi*p2];
H1 = tf(NUM1,DEN1); % transfer function for first pole
H2 = tf(NUM2,DEN2); % transfer function for second pole
H488 = H1 * H2; % transfer function for 2-pole filter
bode(H488); 
[maga,phasea] = bode(H488,w);
for f = 1:Fmax
 Mag488(f) = maga(1,1,f);
 Phase488(f) = phasea(1,1,f);
end
% subplot(2,1,1); 
% plot(Mag488); 
% title('magnitude'); 
% 
% subplot(2,1,2); 
% plot(Phase488); 
% title('phase'); 