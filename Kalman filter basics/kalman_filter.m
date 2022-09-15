clear variables; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates a Kalman filter estimating the position of the hand 
% from simulated noisy data from a single input sensory channel (e.g. vision).
% It also includes an example of sensory fusion where two input channels
% (e.g. vision and proprioception) are combined to yield an estimate of the
% hand's position.
% ----------------------
% Feel free to edit this code and use as you wish!
%
% Atsushi Takagi (2021/03/05) - written and tested in MATLAB 2018b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time step size in seconds
dt = 0.01;

% Length of simulation
Time = 0:dt:2*pi;

% State transition matrix (Euler difference equation)
A = [1 dt dt^2/2;
     0  1     dt;
     0  0     1];

% System noise matrix (this assumes the noise in the position, velocity and
% acceleration are independent of one another). The terms are taken from a
% Taylor expansion of the state transition matrix
Q = diag([dt^3/6, dt^2/2, dt]);

% Observation matrix (this tells you which data you have measurements of).
% In this case, only the position measurement from the hand is available.
C = [1 0 0];

% Noise covariance matrix (this determines the amount of noise in the
% measurements of the hand's position)
Sigma = 2;
R = Sigma^2;

% Initialize state and error covariance matrices
xInit = zeros(3,1);
PInit = diag([1 1 1]).*10^5;

% Generate hand motion and hand's position measurement corrupted by noise
% (e.g. watching the hand move in a dark room)
Signal = sin(Time);
SignalNoise = Signal + sqrt(R)*randn(size(R,1),length(Signal));

% Allocate memory for speed
x = zeros(size(xInit,1),length(Signal));
% Initialize Kalman filter
x(:,1) = xInit;
P = PInit;

% Kalman filter loop
for i=1:length(Signal)-1
    % Prediction step
    x_Prior = A*x(:,i);
    P = A*P*A'+Q;

    % Correction step
    K = P*C'/(C*P*C'+R);
    x(:,i+1) = x_Prior + K*(SignalNoise(:,i)-C*x_Prior);
    P = (eye(size(A,1))-K*C)*P;
end

%% Simulation of sensor fusion between vision and proprioception of hand position

% Sensory measurements come from noisy vision and good proprioception
R2 = diag([Sigma^2 0.05*Sigma^2]);

% Now that two measurements are available, the measurement matrix is
% extended to include the second position measurement from proprioception.
C2 = [1 0 0;
      1 0 0];

% Allocate memory for speed
x2 = zeros(size(xInit,1),length(Signal));
% Initialize Kalman filter
x2(:,1) = xInit;
P = PInit;

% Simulate the noisy measurements of the hand's position
SignalNoise2 = Signal + sqrt(R2)*randn(2,length(Signal));

% Kalman filter loop
for i=1:length(Signal)-1
    % Prediction step
    x_Prior = A*x2(:,i);
    P = A*P*A'+Q;

    % Correction step
    K = P*C2'/(C2*P*C2'+R2);
    x2(:,i+1) = x_Prior + K*(SignalNoise2(:,i)-C2*x_Prior);
    P = (eye(size(A,1))-K*C2)*P;
end

%% PLOTS

figure(1); clf(1); set(gcf,'color','white');

% Figure showing the noisy vision, actual and estimated hand position from
% the Kalman filter. Because vision is very noisy, the estimate is not
% great, albeit it being better than the raw noisy measurements.
subplot(2,1,1); set(gca,'fontsize',15); pbaspect([3,1,1]); hold on;
plot(Time,SignalNoise,'m');
plot(Time,Signal,'k','linewidth',3);
plot(Time,x(1,:),'b','linewidth',2);
legend('noisy vision','actual','estimated','location','northeastoutside'); legend boxoff;
ylabel('position (m)');
ylim([-1,1]*5);

% Figure showing the sensory fusion of noisy vision and good
% proprioception. Because proprioception is good, the estimate improves.
% Note that the estimated value is better than either the noisy vision or
% the good proprioception alone. More sensors and more measurements should
% improve estimation so long as the measurements aren't biased (i.e. offset
% from the actual hand position).
subplot(2,1,2); set(gca,'fontsize',15); pbaspect([3,1,1]); hold on;
plot(Time,SignalNoise2(1,:),'m');
plot(Time,SignalNoise2(2,:),'r','linewidth',1);
plot(Time,Signal,'k','linewidth',3);
plot(Time,x2(1,:),'b','linewidth',2);
legend('noisy vision','good proprio','actual','estimated','location','northeastoutside'); legend boxoff;
xlabel('Time (s)');
ylabel('position (m)');
ylim([-1,1]*5);

