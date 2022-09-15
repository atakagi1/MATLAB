clear variables; clc; close all; set(0,'DefaultFigureWindowStyle','docked'); tic
ThesisPath = 'F:\MEGA\Bioengineering\PhD\PhD Assessments\PhD Thesis\Thesis\figures\svg\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TVINS human-human interaction simulation with different mechanisms
% Connection type
% 1: co-activity
% 2: follow the better
% 3: multi-sensory integration
% 4: interpersonal goal integration
% 5: one-way connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load human-human data
load('interaction_data');

% Which mechanisms to test?
ConnectionArray = 4;

% Number of trials
TotalTrials = size(ImprovementAnalysis,1);
FitTrials = 20;

% Time step size
dt = 0.02;
% Stiffness, damping and mass of partner
k = 120; d = 7; m = 1;

% Visual delay
SensoryDelay = 0.08/dt;

% Generate target in meters
Time=0:dt:30-dt;
Target=(3*sin(1.8*Time)+3.4*sin(1.8*Time)+2.5*sin(1.82*Time)+4.3*sin(2.34*Time))/100;
                   
% Sensory noise parameters (target, hand, partner's hand noise respectively)
MinTaskNoise = 0.005; MaxTaskNoise = 0.1;
TaskNoise = (MinTaskNoise+(MaxTaskNoise-MinTaskNoise)*rand(FitTrials,2)).^2;

% State equation matrices
A = [1  dt dt^2/2;
     0  1  dt;
     0  0  1];
 
% Process noise matrix
Q = 50*[dt^3/6, dt^2/2, dt]'*[dt^3/6, dt^2/2, dt];

% Hand control matrix
B = [0; 0; dt/m];

% Controller state cost
Qdlqr = diag([1,0.1,0]);

% [p_1, v_1, a_1, p_2, v_2, a_2, F, p_t, v_t, a_t, L(1), L(2), L(3)]
InternalModel = @(x)[x(1,:)+x(2,:)*dt+x(3,:)*dt^2/2; x(2,:)+x(3,:)*dt; x(3,:); x(4,:)+x(5,:)*dt+x(6,:)*dt^2/2; x(5,:)+x(6,:)*dt; x(6,:)+((x(7,:)-x(4,:)).*x(11,:)+(x(8,:)-x(5,:)).*x(12,:)+(x(9,:)-x(6,:)).*x(13,:))*B(3); x(7,:)+x(8,:)*dt+x(9,:)*dt^2/2; x(8,:)+x(9,:)*dt; x(9,:); k*(x(4,:)-x(1,:)) + d*(x(5,:)-x(2,:)); x(11,:); x(12,:); x(13,:)];

% Cost of u
Rdlqr = 1*10^-6;
%Rdlqr = 10^-7;

% Compute optimal controller gains
L = dlqr(A,B,Qdlqr,Rdlqr);


ParallelPool = gcp('nocreate');
if isempty(ParallelPool)
    c = parcluster;
    c.NumWorkers = 8;
    parpool('local',8);
end

SingleError = zeros(TotalTrials,2);
DualError = zeros(TotalTrials,2);
CorrTarget = zeros(TotalTrials,2);

% FIT TRIALS TO CONVERT ERROR TO NOISE
parfor Trial=1:FitTrials
    [Robot1,Robot2] = MechanismModel(Target,A,B,Q,L,TaskNoise(Trial,:),SensoryDelay,0,k,d,dt,InternalModel);
    ErrorForFit(Trial,:) = [sqrt(mean((Robot1(1,:)-Target).^2)), sqrt(mean((Robot2(1,:)-Target).^2))];
end
ErrorToTaskNoise = polyfit(ErrorForFit(:),TaskNoise(:),2);
TaskNoiseNew = polyval(ErrorToTaskNoise,ImprovementAnalysis(:,1:2));

% SOLO TRIALS
parfor Trial=1:TotalTrials
    % Interaction simulation (one-way/two-way)
    [Robot1,Robot2] = MechanismModel(Target,A,B,Q,L,TaskNoiseNew(Trial,:),SensoryDelay,0,k,d,dt,InternalModel);
    SingleError(Trial,:) = [sqrt(mean((Robot1(1,:)-Target).^2)), sqrt(mean((Robot2(1,:)-Target).^2))];
end

% CONNECTED TRIALS
for Connection=ConnectionArray
    parfor Trial=1:TotalTrials
        % Interaction simulation (one-way/two-way)
        [Robot1,Robot2] = MechanismModel(Target,A,B,Q,L,TaskNoiseNew(Trial,:),SensoryDelay,Connection,k,d,dt,InternalModel);
        DualError(Trial,:) = [sqrt(mean((Robot1(1,:)-Target).^2)), sqrt(mean((Robot2(1,:)-Target).^2))];
    end
    
    ErrorStore(Connection).DualError = DualError;
end
toc

%% FIGURE: IMPROVEMENT PLOT

ConnectionColor = [1,0,0; 0,0.5,0; 1,0.5,0.1; 0,0,1];
ConnectionString = {'co-activity','follow the better','multi-sensory integration','interpersonal goal integration'};

%%%%%%%%%%% RATIO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scale = 1;

FitOrder = 2;

figure(1); clf(1); set(gcf,'color','w'); set(gca,'fontsize',15); hold on;
line([0 0],[-Scale*2 Scale], 'Color', 'k', 'LineWidth', 0.5);
line([-Scale Scale],[0 0], 'Color', 'k', 'LineWidth', 0.5);

xAxisData = [(ImprovementAnalysis(:,1)-ImprovementAnalysis(:,2))./ImprovementAnalysis(:,1); (ImprovementAnalysis(:,2)-ImprovementAnalysis(:,1))./ImprovementAnalysis(:,2)];
yAxisData = [(ImprovementAnalysis(:,1)-ImprovementAnalysis(:,3))./ImprovementAnalysis(:,1); (ImprovementAnalysis(:,2)-ImprovementAnalysis(:,4))./ImprovementAnalysis(:,2)];

% DATA
[~,~,~,hData,hFill]=HACPlot(xAxisData,yAxisData,'k',FitOrder);
set(hFill,'facealpha',0.1);

% SIMULATIONS
xAxis = [(SingleError(:,1)-SingleError(:,2))./SingleError(:,1); (SingleError(:,2)-SingleError(:,1))./SingleError(:,2)];
for Connection=ConnectionArray
    if Connection==1||Connection==4
        FitOrder = 2;
    else
        FitOrder = 2;
    end
    
    yAxis = (SingleError(:)-ErrorStore(Connection).DualError(:))./SingleError(:);
    [~,~,~,hLine(Connection)]=HACPlot(xAxis,yAxis,ConnectionColor(Connection,:),FitOrder);
    set(hLine(Connection),'linestyle','--');
    
    BetterPartner = yAxis(xAxis<0);
    WorsePartner = yAxis(xAxis>0);
        
end

legend([hData,hLine(ConnectionArray)],cat(2,{'data'},ConnectionString(ConnectionArray)),'interpreter','none','fontsize',8,'location','northeastoutside'); legend boxoff;
axis image; axis([-Scale,Scale/2,-Scale,Scale]); box on;
set(gca,'xtick',[-Scale,0,Scale/2],'ytick',[-Scale,0,Scale]);
xlabel('Partner''s relative error (1-e_p  /e)','fontsize',15);
ylabel('Improvement (1-e_c  /e)','fontsize',15);


