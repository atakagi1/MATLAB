clear variables; close all; clc; set(0,'DefaultFigureWindowStyle','docked');
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program produces the figures used in the following paper:
% ----------------------
% Takagi, A., Li, Y., & Burdet, E. (2020). 
% Flexible Assimilation of Human's Target for Versatile Human-Robot Physical Interaction. 
% IEEE Transactions on Haptics, 14(2), 421-431.
% hand's position.
% ----------------------
% This code showcases the Intention Assimilation Controller (IAC),
% which estimates the human partner's movement intention to enable
% flexible interaction behaviors ranging from collaboration to competition
% using one scalar variable.
%
% Please cite our work if you use this code!
%
% Atsushi Takagi (2023/10/19) - written and tested in MATLAB 2019b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation step size
dt = 0.01;
% Mass
m = 1;
% Duration of simulation
Time = 0:dt:2;

% CHOOSE HUMAN AND ROBOT TARGETS
TargetHuman = 0.3*ones(1,length(Time));

% Target velocity
TargetHuman(2,:) = diff([TargetHuman(1),TargetHuman])/dt;

% System equations
A = [1,dt; 0,1];
B = [0; dt/m];

% Actual controller gain
RHuman = 0.1;
LHuman = dlqr(A,B,diag([100,0*rand(1)]),RHuman)
% Guess
RRobot = 0.1;
LRobot = dlqr(A,B,diag([500,0*rand(1)]),RRobot)
LRobot2 = dlqr(A,B,diag([50,0*rand(1)]),RRobot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up Kalman filter to estimate target from control command u

% Linearization of non-linear state function through complex step differentiation
LinearizeFunc = @(fun,x) imag(fun(x(:,ones(1,numel(x)))+(numel(x)*eps)*1i*eye(numel(x))))/(numel(x)*eps);

AKF = zeros(2,2);
% Observation matrix
H = zeros(1,size(AKF,1));
H(1,end) = 1;

QKF = 10*dt^2;% Process noise

QKF(end+1,end+1) = QKF(end,end);

% Covariance matrix
P = repmat(1000*diag(ones(1,size(AKF,1))),1,1);
P2 = P;

% Noise matrix
RKF = diag(0.001*ones(1,size(H,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulation
uHuman = zeros(1,length(Time));
xHat = zeros(size(AKF,1),length(Time));
xHat2 = xHat;
TargetRobotEstimate = zeros(1,length(Time));
x = zeros(size(A,1),length(Time));
x(:,1) = zeros(size(A,1),1);

% Simulate movement
for i=2:length(Time)-1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TARGET ESTIMATED USING CORRECT GAIN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Internal model to estimate target
    InternalModel = @(xTarget)[xTarget(1,:); (xTarget(1,:)-x(1,i))*LRobot(1,1)+(0-x(2,i))*LRobot(1,2)];
    
    % Linearize internal model
    ATemp = LinearizeFunc(InternalModel,xHat(:,i));

    % Predict
    xTemp=InternalModel(xHat(:,i));
    P=ATemp*P*ATemp'+QKF;
    KTemp=P*H'/(H*P*H'+RKF);
    % Correct
    z = uHuman(:,i-1)+sqrt(RKF)*randn(size(RKF,1),1);
    xHat(:,i+1)=xTemp+KTemp*(z-H*xTemp);
    P=(eye(size(ATemp))-KTemp*H)*P;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TARGET ESTIMATED USING INCORRECT GAIN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Internal model to estimate target
    InternalModel = @(xTarget)[xTarget(1,:); (xTarget(1,:)-x(1,i))*LRobot2(1,1)+(0-x(2,i))*LRobot2(1,2)];
    
    % Linearize internal model
    ATemp = LinearizeFunc(InternalModel,xHat2(:,i));

    % Predict
    xTemp=InternalModel(xHat2(:,i));
    P2=ATemp*P2*ATemp'+QKF;
    KTemp=P2*H'/(H*P2*H'+RKF);
    % Correct
    z = uHuman(:,i-1)+sqrt(RKF)*randn(size(RKF,1),1);
    xHat2(:,i+1)=xTemp+KTemp*(z-H*xTemp);
    P2=(eye(size(ATemp))-KTemp*H)*P2;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Robot control policy
    TargetRobotEstimate(1,i) = xHat(1,i);
    
    % Human control policy
    uHuman(:,i) = -LHuman*(x(:,i)-TargetHuman(:,i));

    % Human and robot system update
    x(:,i+1) = A*x(:,i)+B*uHuman(:,i);
    
end


%% FIGURE: SHOW IDENTIFICATION OF TARGET USING WRONG GAIN
figure(1); clf(1); set(gcf,'color','w'); hold on; set(gca,'fontsize',20,'FontName','Arial'); box on; hold on;

%subplot(2,1,1); set(gca,'fontsize',15,'FontName','Arial'); box on; hold on;
plot(Time,TargetHuman(1,:),'k','linewidth',1);
plot(Time,x(1,:),'k','linewidth',2);
plot(Time,xHat(1,:),'--b','linewidth',2);
plot(Time,xHat2(1,:),'--r','linewidth',2);


set(gca,'xtick',0:0.5:2);
ylim([0,0.5]);
legend('$\tau_h$','$x$','$\hat{\tau}_h\; (||L||<||L_{h}||)$','$\hat{\tau}_h\; (||L||>||L_{h}||)$','location','southeast','interpreter','latex'); legend boxoff;
ylabel('x (m)');
xlabel('Time (s)');
pbaspect([2,1,1]);

%% NEGOTIATION SIMULATION

% CHOOSE HUMAN AND ROBOT TARGETS
TargetHuman = 0.3*ones(1,length(Time));
TargetRobot = -0.3*ones(1,length(Time));

% Actual controller gain
QHuman = diag([100,0*rand(1)]);
LHuman = dlqr(A,B,QHuman,RHuman)
% Guess (once again wrong, but doesn't have to be right as it converges eventually)
QRobot = diag([500,0*rand(1)]);
LRobot = dlqr(A,B,QRobot,RRobot)

% Covariance matrix
PRobot = repmat(1*diag(ones(1,size(AKF,1))),1,1);
PHuman = PRobot;

% Preallocate arrays
xHatHuman = xHat;
xHatRobot = xHat;
xCoactivity = x;
xHumanSolo = x;
xRobotSolo = x;
xPassiveFollow = x;
xNegotiate = zeros(size(A,1),length(Time),3);
uRobotNegotiate = zeros(1,length(Time),3); uHumanNegotiate = uRobotNegotiate;
uRobot = uHuman;

% Simulate movement for different scenarios
for s=1:3
    
    switch s
        case 1
            % Robot assists human to reach their target
            LambdaRobot = 0;
            LambdaHuman = 1;
        case 2
            % Coactivity human and robot each their targets, essentially
            % ignoring each other. Final position will depend on their
            % relative controller gain size.
            LambdaRobot = 1;
            LambdaHuman = 1;
        case 3
            % Robot competes with human to reach the robot's target
            LambdaRobot = 2;
            LambdaHuman = 1;
    end
    
    for i=2:length(Time)-1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TARGET ESTIMATED USING ROBOT CONTROLLER GAIN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Internal model to estimate target
        InternalModel = @(xTarget)[xTarget(1,:); (xTarget(1,:)-xNegotiate(1,i,s))*LHuman(1,1)+(0-xNegotiate(2,i,s))*LHuman(1,2)];

        % Linearize internal model
        ATemp = LinearizeFunc(InternalModel,xHatHuman(:,i));

        % Predict
        xTemp=InternalModel(xHatHuman(:,i));
        PHuman=ATemp*PHuman*ATemp'+QKF;
        KTemp=PHuman*H'/(H*PHuman*H'+RKF);
        % Correct
        z = uRobotNegotiate(:,i-1,s)+sqrt(RKF)*randn(size(RKF,1),1);
        xHatHuman(:,i+1)=xTemp+KTemp*(z-H*xTemp);
        PHuman=(eye(size(ATemp))-KTemp*H)*PHuman;
        %%%%%%%%%%%%%%%%%%%%%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TARGET ESTIMATED USING ROBOT CONTROLLER GAIN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Internal model to estimate target
        InternalModel = @(xTarget)[xTarget(1,:); (xTarget(1,:)-xNegotiate(1,i,s))*LRobot(1,1)+(0-xNegotiate(2,i,s))*LRobot(1,2)];

        % Linearize internal model
        ATemp = LinearizeFunc(InternalModel,xHatRobot(:,i));

        % Predict
        xTemp=InternalModel(xHatRobot(:,i));
        P=ATemp*P*ATemp'+QKF;
        KTemp=P*H'/(H*P*H'+RKF);
        % Correct
        z = uHumanNegotiate(:,i-1,s)+sqrt(RKF)*randn(size(RKF,1),1);
        xHatRobot(:,i+1)=xTemp+KTemp*(z-H*xTemp);
        P=(eye(size(ATemp))-KTemp*H)*P;
        %%%%%%%%%%%%%%%%%%%%%%%    

        % Target update
        TargetRobotEstimate(1,i) = (LambdaRobot*TargetRobot(1,i)+(1-LambdaRobot)*xHatRobot(1,i));
        TargetHumanEstimate(1,i) = (LambdaHuman*TargetHuman(1,i)+(1-LambdaHuman)*xHatHuman(1,i));


        % Control policy
        uRobotNegotiate(:,i,s) = -LRobot*(xNegotiate(:,i,s)-[TargetRobotEstimate(1,i); 0]);
        uHumanNegotiate(:,i,s) = -LHuman*(xNegotiate(:,i,s)-[TargetHumanEstimate(1,i); 0]);

        % Human and robot system update
        xNegotiate(:,i+1,s) = A*xNegotiate(:,i,s)+B*(uHumanNegotiate(:,i,s)+uRobotNegotiate(:,i,s));
    end
end
    
%% FIGURE: NEGOTIATION
figure(2); clf(2); set(gcf,'color','w'); hold on; set(gca,'fontsize',20,'FontName','Arial'); box on; hold on;

%subplot(2,1,1); set(gca,'fontsize',15,'FontName','Arial'); box on; hold on;
plot(Time,TargetHuman(1,:),'k','linewidth',1);
plot(Time,TargetRobot(1,:),'--k','linewidth',1);
h(1)=plot(Time,xNegotiate(1,:,1),'color',[0,0.5,0],'linewidth',2);
h(2)=plot(Time,xNegotiate(1,:,2),'k','linewidth',2); % coactivity
h(3)=plot(Time,xNegotiate(1,:,3),'r','linewidth',2);

text(0.2,TargetHuman(1,1)*1.12,'$\tau_h$','interpreter','latex','fontsize',20);
text(0.2,TargetRobot(1,1)*1.12,'$\tau_r$','interpreter','latex','fontsize',20);

set(gca,'xtick',0:0.5:2);
ylim([-0.4,0.4]);
hLegend=legend(h,'Assistance ($\lambda<1$)','Coactivity ($\lambda=1$)','Antagonistic behaviour ($\lambda>1$)','location','northoutside','interpreter','latex'); legend boxoff;
rect = [0.55, 0.47, .1, .17];
set(hLegend, 'Position', rect);
ylabel('x (m)');
xlabel('Time (s)');
pbaspect([2,1,1]);

%% SIMULATE HUMAN REACHING ALONE, OR WITH COLLABORATIVE IAC, OR WITH PASSIVE FOLLOWER

for i=2:length(Time)-1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    uRobot = -LRobot*(xPassiveFollow(:,i)-[xPassiveFollow(1,i); 0]);
    uHumanPassiveFollow(1,i) = -LHuman*(xPassiveFollow(:,i)-[TargetHuman(1,i); 0]);
    xPassiveFollow(:,i+1) = A*xPassiveFollow(:,i)+B*(uRobot+uHumanPassiveFollow(1,i));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Human alone system update
    uHumanSolo(:,i) = -LHuman*(xHumanSolo(:,i)-[TargetHuman(1,i);0]);
    xHumanSolo(:,i+1) = A*xHumanSolo(:,i)+B*uHumanSolo(:,i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Robot alone system update
    uRobotSolo = -LRobot*(xRobotSolo(:,i)-[TargetRobot(1,i);0]);
    xRobotSolo(:,i+1) = A*xRobotSolo(:,i)+B*uRobotSolo;
    
end
%% FIGURE: COMPARISON WITH OTHER STRATEGIES

figure(3); clf(3); set(gcf,'color','w'); hold on; set(gca,'fontsize',15,'FontName','Arial'); box on; hold on;

subplot(2,1,1); set(gca,'fontsize',20,'FontName','Arial'); box on; hold on;

%plot(Time,TargetRobot(1,:),'--m','linewidth',1);
plot(Time,xHumanSolo(1,:),'b','linewidth',2);
plot(Time,xNegotiate(1,:,1),'color',[0,0.5,0],'linewidth',2);
plot(Time,xPassiveFollow(1,:),'k','linewidth',2);
plot(Time,TargetHuman(1,:),'k','linewidth',1);

text(0.2,TargetHuman(1,1)*1.06,'$\tau_h$','interpreter','latex','fontsize',20);
set(gca,'xtick',0:0.5:2);
ylim([0,0.4]);
legend('human solo','$\lambda=0$','passive follower','location','best','interpreter','latex'); legend boxoff;
ylabel('x (m)');
xlabel('Time (s)');


subplot(2,1,2); set(gca,'fontsize',20,'FontName','Arial'); box on; hold on;
plot(Time(2:end-1),uHumanSolo(1,2:end),'b','linewidth',2);
plot(Time(2:end-2),uHumanNegotiate(1,2:end-2,1),'color',[0,0.5,0],'linewidth',2);
plot(Time(2:end-1),uHumanPassiveFollow(1,2:end),'k','linewidth',2);
line([0,2],[0,0],'color','k');
legend('human solo','$\lambda=0$','passive follower','location','best','interpreter','latex'); legend boxoff;
set(gca,'xtick',0:0.5:2);
ylim([-5,10]);
ylabel('u_h (N)');
xlabel('Time (s)');

%% COLLISION AVOIDANCE VIA COMPETITIVE IAC

xWall = 0.06;

dt = 0.01;

Time = 0:dt:3;

% System equations
A = [1,dt; 0,1];
B = [0; dt/m];


QHuman = diag([150,0]);
LHuman = dlqr(A,B,QHuman,RHuman)
% Guess
QRobot = diag([500,0]);
LRobot = dlqr(A,B,QRobot,RRobot)

QRobot = diag([80,0]);
LRobot2 = dlqr(A,B,QRobot,RRobot)


xCollision = zeros(2,length(Time));
xCollision2 = xCollision;
xCollisionHomotopy = xCollision;
LambdaRobotCollision = zeros(1,length(Time)-1);
LambdaRobotCollision2 = LambdaRobotCollision;
uHumanCollision = zeros(1,length(Time)-1);
uRobotCollision = uHumanCollision;

% Process noise
QKF = 10*dt^2;
QKF(end+1,end+1) = QKF(end,end);

RKF = 0.00001;
xHatRobot = xHat;
xHatRobot2 = xHat;
P2 = P;
for i=2:length(Time)-1
    
    % This is human's target
    if i<1/dt
        TargetHumanCollide(1,i) = (0.15/1) * i*dt;
    else
        TargetHumanCollide(1,i) = max(0,-(0.15/1) * (i*dt-1) + 0.15);
    end       

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TARGET ESTIMATED USING ROBOT CONTROLLER GAIN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Internal model to estimate target
    InternalModel = @(xTarget)[xTarget(1,:); (xTarget(1,:)-xCollision(1,i))*LRobot(1,1)+(0-xCollision(2,i))*LRobot(1,2)];

    % Linearize internal model
    ATemp = LinearizeFunc(InternalModel,xHatRobot(:,i));

    % Predict
    xTemp=InternalModel(xHatRobot(:,i));
    P=ATemp*P*ATemp'+QKF;
    KTemp=P*H'/(H*P*H'+RKF);
    % Correct
    z = uHumanCollision(:,i-1)+sqrt(RKF)*randn(size(RKF,1),1);
    xHatRobot(:,i+1)=xTemp+KTemp*(z-H*xTemp);
    P=(eye(size(ATemp))-KTemp*H)*P;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Internal model to estimate target
    InternalModel = @(xTarget)[xTarget(1,:); (xTarget(1,:)-xCollision2(1,i))*LRobot2(1,1)+(0-xCollision2(2,i))*LRobot2(1,2)];

    % Linearize internal model
    ATemp = LinearizeFunc(InternalModel,xHatRobot2(:,i));

    % Predict
    xTemp=InternalModel(xHatRobot2(:,i));
    P2=ATemp*P2*ATemp'+QKF;
    KTemp=P2*H'/(H*P2*H'+RKF);
    % Correct
    z = uHumanCollision(:,i-1)+sqrt(RKF)*randn(size(RKF,1),1);
    xHatRobot2(:,i+1)=xTemp+KTemp*(z-H*xTemp);
    P2=(eye(size(ATemp))-KTemp*H)*P2;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    LambdaRobotCollision(1,i) = max(0,min(2,2/(0.05)*(xHatRobot(1,i)+0.05-xWall)));
    LambdaRobotCollision2(1,i) = max(0,min(2,2/(0.05)*(xHatRobot2(1,i)+0.05-xWall)));
    AlphaRobotCollision(1,i) = max(0,min(2,2/(0.05)*(xCollisionHomotopy(1,i)+0.05-xWall)));

    % Target update
    TargetRobotEstimate(1,i) = (LambdaRobotCollision(1,i)*(xCollision(1,i))+(1-LambdaRobotCollision(1,i))*xHatRobot(1,i));
    TargetRobotEstimate2(1,i) = (LambdaRobotCollision2(1,i)*(xCollision2(1,i))+(1-LambdaRobotCollision2(1,i))*xHatRobot2(1,i));

    % Control policy
    uRobotCollision(:,i) = -LRobot*(xCollision(:,i)-[TargetRobotEstimate(1,i); 0]);
    uHumanCollision(:,i) = -LHuman*(xCollision(:,i)-[TargetHumanCollide(1,i); 0]);
    uRobotCollision2(:,i) = -LRobot2*(xCollision2(:,i)-[TargetRobotEstimate2(1,i); 0]);
    uHumanCollision2(:,i) = -LHuman*(xCollision2(:,i)-[TargetHumanCollide(1,i); 0]);
    uHumanCollisionHomotopy(:,i) = -LHuman*(xCollisionHomotopy(:,i)-[TargetHumanCollide(1,i); 0]);
    uRobotCollisionHomotopy(:,i) = AlphaRobotCollision(1,i)*-LRobot*(xCollisionHomotopy(:,i)-[xCollisionHomotopy(1,i); 0]) + (1-AlphaRobotCollision(1,i))*uHumanCollisionHomotopy(:,i);

    % Human and robot system update
    xCollision(:,i+1) = A*xCollision(:,i)+B*(uHumanCollision(:,i)+uRobotCollision(:,i));
    xCollision2(:,i+1) = A*xCollision2(:,i)+B*(uHumanCollision2(:,i)+uRobotCollision2(:,i));
    xCollisionHomotopy(:,i+1) = A*xCollisionHomotopy(:,i)+B*(uRobotCollisionHomotopy(:,i) + uHumanCollisionHomotopy(:,i));
end

%% FIGURE: HOMOTOPY DOES NOT PREVENT COLLISION WITH WALL, BUT COMPETITIVE IAC DOES

figure(4); clf(4); set(gcf,'color','w'); hold on; set(gca,'fontsize',15,'FontName','Arial'); box on; hold on;

subplot(2,1,1); set(gca,'fontsize',20,'FontName','Arial'); box on; hold on;
plot(Time,xWall*ones(1,length(Time)),'r','linewidth',2);
hCollision(1)=plot(Time,xHatRobot(1,:),'color','m','linestyle','-','linewidth',2);
hCollision(2)=plot(Time,xHatRobot2(1,:),'color','m','linestyle',':','linewidth',2);
hCollision(3)=plot(Time,xCollision(1,:),'k','linewidth',2);
hCollision(4)=plot(Time,xCollision2(1,:),':k','linewidth',2);
hCollision(5)=plot(Time,xCollisionHomotopy(1,:),'b','linewidth',2);

plot(Time(1:end-1),TargetHumanCollide,'k','linewidth',0.5);

text(0.11,xWall*1.4,'$\tau_{\textup{wall}}$','interpreter','latex','fontsize',20,'color','r');
hText=text(0.8,0.11,'$\tau_h$','interpreter','latex','fontsize',20);
set(hText,'Rotation',25);
set(gca,'xtick',0:0.5:3);
ylim([-0.02,0.25]);
legend(hCollision,'$\hat{\tau}_{h}\,(||L^{v}_h||>||L_h||)$','$\hat{\tau}_{h}\,(||L^{v}_h||<||L_h||)$','$\textup{IAC}\,(||L^{v}_h||>||L_h||)$','$\textup{IAC}\,(||L^{v}_h||<||L_h||)$','$\textup{Homotopy}$','interpreter','latex'); legend boxoff;
ylabel('x (m)');
subplot(2,1,2); set(gca,'fontsize',20,'FontName','Arial'); box on; hold on;
plot(Time(2:end-1),LambdaRobotCollision(1,2:end),'k','linewidth',2);
plot(Time(2:end-1),LambdaRobotCollision2(1,2:end),':k','linewidth',2);
legend('$||L^{v}_h||>||L_h||$','$||L^{v}_h||<||L_h||$','interpreter','latex'); legend boxoff;
set(gca,'xtick',0:0.5:3);
ylim([-0.5,2.5]);
ylabel({'$\lambda$'},'Interpreter','latex');
xlabel('Time (s)');


