clear variables; clc; close all; set(0,'DefaultFigureWindowStyle','docked'); tic
set(0,'DefaultAxesFontName','Arial');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of two partners tracking a randomly moving target with a
% rigid, medium or compliant physical coupling
% 
% This simulation is split into two parts:
% 1. Simulate solo trials (no interaction) to obtain the relationship
% between error and noise for the system (this relationship depends on the
% control gain and the system noise matrix in the Kalman filter). This fit
% is used to set the sensory noise needed to produce partners within a
% specified error range, and it is necessary to calculate the additional
% noise that comes from the compliance in the physical coupling. When the
% coupling is rigid, the partner's goal has little additional noise, but
% when the coupling is compliant it is corrupted by significant noise. 
%
% 2. Simulate interaction trials. Three seperate instances are conducted
% (for rigid, medium and compliant coupling). The errors and the effort are
% all stored for later analysis.
%
% The simulations are compared with the data stored in the 'HHstiffdata.mat'
% file.
%
% Atsushi Takagi (last updated 2023/04/11)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of trials to simulate
TotalTrials = 20;

% Connection type (1=relaxation,2=two-way,3=naive two-way,4=coactivity,5=follow the better)
ConnectionArray = 2;

% Time step size
dt = 0.02;
% Inertia [kg m^2] (set heavy to match mean RA torque in humans)
I = 0.015;
% Cost of u
Rdlqr = 10^-10;

% Generate target in metres
Time=0:dt:40-dt;
Target = 1.7*(-2.63*(-2.376*sin(-3.061*Time*1.3)+0.366*sin(1.347*Time*1.3)-1.734*sin(-0.579*Time*1.3)-2.092*sin(1.774*Time*1.3)));
                   
% Task noise [deg^2] (converted to rad^2 in loop) (uniformly distribute noise STD, not noise)
MinTaskNoiseSTD = 0.1; MaxTaskNoiseSTD = 5;
TaskNoiseVisual = (MinTaskNoiseSTD+(MaxTaskNoiseSTD-MinTaskNoiseSTD)*rand(TotalTrials,2)).^2;

% System equation
A = [1, dt;
     0, 1];

% Hand control matrix
B = [0; dt/I];

% Allocate memory
SingleError1 = zeros(TotalTrials,3);
SingleError2 = zeros(TotalTrials,3);
EffortSingle1 = zeros(TotalTrials,3);
EffortSingle2 = zeros(TotalTrials,3);

% Control gain modifier (used in sensitivity analysis)
q = [0.001,0.01];
% Base system noise matrix
QBase = [dt^2/2, dt]'*[dt^2/2, dt];
% System noise matrix modifier (used in sensitivity analysis)
mu = [1,10,50];

SimData = cell(length(q),length(mu),max(ConnectionArray));
MetricError = zeros(length(q),length(mu),max(ConnectionArray)); MetricEffort = MetricError;
for m=1:length(q)

    % Controller state cost
    Qdlqr = diag([1,q(m)]);
    % Compute optimal controller gains
    L = dlqr(A,B,Qdlqr,Rdlqr)

    for n=1:length(mu)

        % System noise matrix
        Q = mu(n)*QBase;
        
        % Simulate solo trials to get relationship between error and noise
        parfor Trial=1:TotalTrials
            % Solo/single simulation
            [SingleRobot1,SingleRobot2,u1,u2] = compliance_model2(Target*pi/180,A,B,Q,L,TaskNoiseVisual(Trial,:)*(pi/180)^2,0,0,dt,I);
            % Solo/single error
            SingleError(Trial,:) = [sqrt(mean((SingleRobot1(1,:)-Target(1,:)).^2)), sqrt(mean((SingleRobot2(1,:)-Target(1,:)).^2))];
        
        end
        % Store
        SimData{m,n,1}.SingleError = SingleError(:);
        SimData{m,n,1}.TaskNoise = TaskNoiseVisual(:);
        
        % Fit error and noise as functions of each other
        ErrorToNoiseParam = polyfit(SingleError(:),sqrt(TaskNoiseVisual(:)),2);
        NoiseToErrorParam = polyfit(sqrt(TaskNoiseVisual(:)),SingleError(:),2);

        SimData{m,n,1}.ErrorToNoiseParam = ErrorToNoiseParam;

        % Error range for partners
        MinError = 5; MaxError = 12;
        TaskNoise = max(0.001,polyval(ErrorToNoiseParam,(MinError+(MaxError-MinError)*rand(TotalTrials,2))).^2);

        for Trial=1:TotalTrials
            % Solo/single simulation
            [SingleRobot1,SingleRobot2,u1,u2] = compliance_model2(Target*pi/180,A,B,Q,L,TaskNoise(Trial,:)*(pi/180)^2,0,0,dt,I);
            % Solo/single error
            SingleError1(Trial,:) = sqrt(mean((SingleRobot1(1,:)-Target(1,:)).^2))*ones(1,3); SingleError2(Trial,:) = sqrt(mean((SingleRobot2(1,:)-Target(1,:)).^2))*ones(1,3);
            EffortSingle1(Trial,:) = mean(abs(u1))*ones(1,3); EffortSingle2(Trial,:) = mean(abs(u2))*ones(1,3);
        end

        DualErrorRigid1=zeros(TotalTrials,1); DualErrorRigid2=DualErrorRigid1;
        DualErrorMed1=DualErrorRigid1; DualErrorMed2=DualErrorRigid1;
        DualErrorLow1=DualErrorRigid1; DualErrorLow2=DualErrorRigid1;

        DualEffortRigid1=zeros(TotalTrials,1); DualEffortRigid2=DualEffortRigid1;
        DualEffortMed1=DualEffortRigid1; DualEffortMed2=DualEffortRigid1;
        DualEffortLow1=DualEffortRigid1; DualEffortLow2=DualEffortRigid1;

        for Connection=ConnectionArray
            for Trial=1:TotalTrials

                % Interaction simulation (rigid)
                [DualRobot1,DualRobot2,u1,u2] = compliance_model2(Target*pi/180,A,B,Q,L,TaskNoise(Trial,:)*(pi/180)^2,Connection,0.3*180/pi,ErrorToNoiseParam,NoiseToErrorParam);
                DualErrorRigid1(Trial,1) = sqrt(mean((DualRobot1(1,:)-Target(1,:)).^2)); DualErrorRigid2(Trial,1) = sqrt(mean((DualRobot2(1,:)-Target(1,:)).^2));
                DualEffortRigid1(Trial,1) = mean(abs(u1)); DualEffortRigid2(Trial,1) = mean(abs(u2));

                % Interaction simulation (medium)
                [DualRobot1,DualRobot2,u1,u2] = compliance_model2(Target*pi/180,A,B,Q,L,TaskNoise(Trial,:)*(pi/180)^2,Connection,0.03*180/pi,ErrorToNoiseParam,NoiseToErrorParam);
                DualErrorMed1(Trial,1) = sqrt(mean((DualRobot1(1,:)-Target(1,:)).^2)); DualErrorMed2(Trial,1) = sqrt(mean((DualRobot2(1,:)-Target(1,:)).^2));
                DualEffortMed1(Trial,1) = mean(abs(u1)); DualEffortMed2(Trial,1) = mean(abs(u2));

                % Interaction simulation (compliant)
                [DualRobot1,DualRobot2,u1,u2] = compliance_model2(Target*pi/180,A,B,Q,L,TaskNoise(Trial,:)*(pi/180)^2,Connection,0.005*180/pi,ErrorToNoiseParam,NoiseToErrorParam);
                DualErrorLow1(Trial,1) = sqrt(mean((DualRobot1(1,:)-Target(1,:)).^2)); DualErrorLow2(Trial,1) = sqrt(mean((DualRobot2(1,:)-Target(1,:)).^2));
                DualEffortLow1(Trial,1) = mean(abs(u1)); DualEffortLow2(Trial,1) = mean(abs(u2));

            end

            % Store all simulation results
            SimData{m,n,Connection}.SingleError1 = SingleError1;
            SimData{m,n,Connection}.SingleError2 = SingleError2;
            SimData{m,n,Connection}.SingleEffort1 = EffortSingle1;
            SimData{m,n,Connection}.SingleEffort2 = EffortSingle2;
            SimData{m,n,Connection}.DualError1 = [DualErrorRigid1,DualErrorMed1,DualErrorLow1];
            SimData{m,n,Connection}.DualError2 = [DualErrorRigid2,DualErrorMed2,DualErrorLow2];
            SimData{m,n,Connection}.DualEffort1 = [DualEffortRigid1,DualEffortMed1,DualEffortLow1];
            SimData{m,n,Connection}.DualEffort2 = [DualEffortRigid2,DualEffortMed2,DualEffortLow2];
        end


        

    end
end

SimData{1,1,1}.q = q;
SimData{1,1,1}.mu = mu;

toc


%%
% PARAMETERS
Color = [1,0,0; 0,0,1; 0,0.6,0];

load('HHstiffdata.mat');
% COLUMN ID FROM 'HHstiffdata'
ColSErrorG = 1;
ColSErrorY = 2;
ColDErrorG = 3;
ColDErrorY = 4;
ColSEffortG = 5;
ColSEffortY = 6;
ColDEffortG = 7;
ColDEffortY = 8;
ColSCoconG = 9;
ColSCoconY = 10;
ColDCoconG = 11;
ColDCoconY = 12;
ColSRecipG = 13;
ColSRecipY = 14;
ColDRecipG = 15;
ColDRecipY = 16;
ColStiffness = 17;
ColDyad = 18;
ColSex = 20; % 0: male, 1: female

% Get data from connection stiffness
SelectStiff = @(Data,ColData,Stiffness)Data(Data(:,ColStiffness)==Stiffness,ColData);

StiffArray = [0.3,0.03,0.005];

SimLMELabel = {'Diff','Diff2','Improve','Coupling'};


x = transpose(-1:0.01:0.5);

%% PERFORMANCE IMPROVEMENT CALCULATION

ScaleError = 1;

SimErrorLME = cell(length(SimData{1,1,1}.q),length(SimData{1,1,1}.mu),max(ConnectionArray));
SimErrorPoly = cell(length(SimData{1,1,1}.q),length(SimData{1,1,1}.mu),max(ConnectionArray),3);

% Calculate points from data to compare with simulation (only needed once)
for Stiff=1:3

    xData{Stiff} = [(SelectStiff(DataStore,ColSErrorG,StiffArray(Stiff))-SelectStiff(DataStore,ColSErrorY,StiffArray(Stiff)))./SelectStiff(DataStore,ColSErrorG,StiffArray(Stiff)); (SelectStiff(DataStore,ColSErrorY,StiffArray(Stiff))-SelectStiff(DataStore,ColSErrorG,StiffArray(Stiff)))./SelectStiff(DataStore,ColSErrorY,StiffArray(Stiff))];
    yErrorData{Stiff} = [(SelectStiff(DataStore,ColSErrorG,StiffArray(Stiff))-SelectStiff(DataStore,ColDErrorG,StiffArray(Stiff)))./SelectStiff(DataStore,ColSErrorG,StiffArray(Stiff)); (SelectStiff(DataStore,ColSErrorY,StiffArray(Stiff))-SelectStiff(DataStore,ColDErrorY,StiffArray(Stiff)))./SelectStiff(DataStore,ColSErrorY,StiffArray(Stiff))];

    yEffortData{Stiff} = (cat(1,(SelectStiff(DataStore,ColDEffortG,StiffArray(Stiff))-SelectStiff(DataStore,ColSEffortG,StiffArray(Stiff)))./SelectStiff(DataStore,ColSEffortG,StiffArray(Stiff)),(SelectStiff(DataStore,ColDEffortY,StiffArray(Stiff))-SelectStiff(DataStore,ColSEffortY,StiffArray(Stiff)))./SelectStiff(DataStore,ColSEffortY,StiffArray(Stiff))));
    
    [yTemp,yTempSigma]=predict(HHLME2,table(xData{Stiff},xData{Stiff}.^2,zeros(length(xData{Stiff}),1),zeros(length(xData{Stiff}),1),log10(StiffArray(Stiff)*ones(length(xData{Stiff}),1)),'VariableNames',{'Diff','Diff2','Improve','Dyad','Coupling'}));
    
    y_data_fit_error{Stiff} = yTemp;
    y_data_var_error{Stiff} = (yTemp-yTempSigma(:,1)).^2;
    
    [yTemp,yTempSigma]=predict(HHEffortLME2,table(xData{Stiff},xData{Stiff}.^2,zeros(length(xData{Stiff}),1),zeros(length(xData{Stiff}),1),log10(StiffArray(Stiff)*ones(length(xData{Stiff}),1)),'VariableNames',LMELabel));
    
    y_data_fit_effort{Stiff} = yTemp;
    y_data_var_effort{Stiff} = (yTemp-yTempSigma(:,1)).^2;
end

for Connection=ConnectionArray
    for m=1:length(SimData{1,1,1}.q)
        for n=1:length(SimData{1,1,1}.mu)
            

            SimDiff = [(SimData{m,n,Connection}.SingleError1(:)-SimData{m,n,Connection}.SingleError2(:))./SimData{m,n,Connection}.SingleError1(:); (SimData{m,n,Connection}.SingleError2(:)-SimData{m,n,Connection}.SingleError1(:))./SimData{m,n,Connection}.SingleError2(:)];
            SimImprove = [(SimData{m,n,Connection}.SingleError1(:)-SimData{m,n,Connection}.DualError1(:))./SimData{m,n,Connection}.SingleError1(:); (SimData{m,n,Connection}.SingleError2(:)-SimData{m,n,Connection}.DualError2(:))./SimData{m,n,Connection}.SingleError2(:)];
            Coupling = repmat(log10([0.3*ones(size(SimData{m,n,Connection}.SingleError1,1),1); 0.03*ones(size(SimData{m,n,Connection}.SingleError1,1),1); 0.005*ones(size(SimData{m,n,Connection}.SingleError1,1),1)]),2,1);
            CouplingNonlog = repmat(([0.3*ones(size(SimData{m,n,Connection}.SingleError1,1),1); 0.03*ones(size(SimData{m,n,Connection}.SingleError1,1),1); 0.005*ones(size(SimData{m,n,Connection}.SingleError1,1),1)]),2,1);
 
            % FOR CALCULATING
            SimErrorLME{m,n,Connection} = fitlme(dataset({cat(2,SimDiff,SimDiff.^2,SimImprove,Coupling),SimLMELabel{:}}),'Improve ~ 1 + Diff*Coupling + Diff2*Coupling','FitMethod','ML');
            
            for Stiff=1:3

                % Simulation
               [y_sim_fit,~]=predict(SimErrorLME{m,n,Connection},dataset({[xData{Stiff},xData{Stiff}.^2,zeros(length(xData{Stiff}),1),log10(StiffArray(Stiff)*ones(length(xData{Stiff}),1))],SimLMELabel{:}}));
                
                % Calculate SSE for error
                SimData{m,n,Connection}.MetricError(1,Stiff) = mean(abs(y_sim_fit-y_data_fit_error{Stiff})./sqrt(y_data_var_error{Stiff}));
            end

            MetricError(m,n,Connection) = nanmean(SimData{m,n,Connection}.MetricError);

        end
    end
end


%% VISUAL NOISE VERSUS TRACKING ERROR

figure(2); close(2); f2=figure(2); set(gcf,'color','w'); set(gca,'fontsize',15);

xError = 1:0.1:10;

for m=1:length(SimData{1,1,1}.q)
    for n=1:length(SimData{1,1,1}.mu)
        
        subplot(length(SimData{1,1,1}.q),length(SimData{1,1,1}.mu),(m-1)*length(SimData{1,1,1}.mu)+n); hold on;
        
        plot(SimData{m,n,1}.SingleError,sqrt(SimData{m,n,1}.TaskNoise),'.b');
        plot(xError,polyval(SimData{m,n,1}.ErrorToNoiseParam,xError),'b');
        
        xlabel('error');
        ylabel('noise');
    end
end

%% EFFORT INCREASE CALCULATION

if Connection==1
    ScaleEffort = 10;
elseif Connection==2||Connection==3
    ScaleEffort = 1;
elseif Connection==4
    ScaleEffort = 2;
elseif Connection==5
    ScaleEffort = 4;
end

SimEffortLME = cell(length(SimData{1,1,1}.q),length(SimData{1,1,1}.mu),max(ConnectionArray));
SimEffortPoly = cell(length(SimData{1,1,1}.q),length(SimData{1,1,1}.mu),max(ConnectionArray),3);

SimEffortLMELabel = {'Diff','Diff2','Diff3','Improve','Coupling'};

for Connection=ConnectionArray
    for m=1:length(SimData{1,1,1}.q)
        for n=1:length(SimData{1,1,1}.mu)

            SimDiff = [(SimData{m,n,Connection}.SingleError1(:)-SimData{m,n,Connection}.SingleError2(:))./SimData{m,n,Connection}.SingleError1(:); (SimData{m,n,Connection}.SingleError2(:)-SimData{m,n,Connection}.SingleError1(:))./SimData{m,n,Connection}.SingleError2(:)];
            SimImprove = [(SimData{m,n,Connection}.SingleError1(:)-SimData{m,n,Connection}.DualError1(:))./SimData{m,n,Connection}.SingleError1(:); (SimData{m,n,Connection}.SingleError2(:)-SimData{m,n,Connection}.DualError2(:))./SimData{m,n,Connection}.SingleError2(:)];
            SimEffort = [(SimData{m,n,Connection}.DualEffort1(:)-SimData{m,n,Connection}.SingleEffort1(:))./SimData{m,n,Connection}.SingleEffort1(:); (SimData{m,n,Connection}.DualEffort2(:)-SimData{m,n,Connection}.SingleEffort2(:))./SimData{m,n,Connection}.SingleEffort2(:)];
            
            SimEffortLME{m,n,Connection} = fitlme(dataset({cat(2,SimDiff,SimDiff.^2,SimDiff.^3,SimEffort,Coupling),SimEffortLMELabel{:}}),'Improve ~ 1 + Diff*Coupling','FitMethod','ML');

            
            for Stiff=1:3
                
                % Simulation
                [y_sim_fit,~]=predict(SimEffortLME{m,n,Connection},dataset({[xData{Stiff},xData{Stiff}.^2,xData{Stiff}.^3,zeros(length(xData{Stiff}),1),log10(StiffArray(Stiff)*ones(length(xData{Stiff}),1))],SimEffortLMELabel{:}}));

                
                % Calculate SSE for effort               
                SimData{m,n,Connection}.MetricEffort(1,Stiff) = mean(abs(y_sim_fit-y_data_fit_effort{Stiff})./sqrt(y_data_var_effort{Stiff}));
            end

            MetricEffort(m,n,Connection) = nanmean(SimData{m,n,Connection}.MetricEffort);
            

        end
    end
end


%% FIGURE: METRIC PLOT FOR ALL STRATEGIES

MetricData = 0.5*MetricError+0.5*MetricEffort;

f=figure(3); clf(3); set(gcf,'color','none'); set(gca,'fontsize',15); hold on;

RescaleMetric = @(Metric)log10(Metric);

for Connection=ConnectionArray
    MetricData(:,:,Connection) = RescaleMetric(MetricData(:,:,Connection));
    MetricData(~isfinite(MetricData(:,:,Connection)),Connection) = NaN;
    
    surf(SimData{1,1,1}.mu,SimData{1,1,1}.q,MetricData(:,:,Connection));
    xlabel('\sigma^2_{\mu}');
    ylabel('q');    
    zlabel('log_{10}(MAE)');
end

view(-10,10); grid on; box on;

ylim([min(SimData{1,1,1}.q),max(SimData{1,1,1}.q)]);
xlim([min(SimData{1,1,1}.mu),max(SimData{1,1,1}.mu)]);

%% FIGURE: IMPROVEMENT + EFFORT FROM LOWEST METRIC FOR ALL STRATEGIES

SubplotFontSize = 10;

for Connection=ConnectionArray
    
    if Connection==1
        ScaleEffort = 10;
    elseif Connection==2||Connection==3
        ScaleEffort = 1;
    elseif Connection==4
        ScaleEffort = 2;
    elseif Connection==5
        ScaleEffort = 6;
    end
    
    figure(50+Connection); close(50+Connection); f=figure(50+Connection); set(gcf,'color','w'); hold on;
    
    TempMetricData = MetricData(:,:,Connection);
    MinMetricData = min(TempMetricData(:));
    [MinM,MinN] = find(MetricData(:,:,Connection)==MinMetricData);

    m = MinM; n = MinN;

    for Stiff=1:3

        subplot(1,2,1); set(gca,'fontsize',SubplotFontSize); hold on;
        line([-ScaleError ScaleError],[0 0], 'Color', 'k', 'LineWidth', 1);
        line([0 0],[-ScaleError*10 ScaleError], 'Color', 'k', 'LineWidth', 1);
        
        % Simulation plot
        %plot([(SimData{m,n,Connection}.SingleError1(:,Stiff)-SimData{m,n,Connection}.SingleError2(:,Stiff))./SimData{m,n,Connection}.SingleError1(:,Stiff); (SimData{m,n,Connection}.SingleError2(:,Stiff)-SimData{m,n,Connection}.SingleError1(:,Stiff))./SimData{m,n,Connection}.SingleError2(:,Stiff)],[(SimData{m,n,Connection}.SingleError1(:,Stiff)-SimData{m,n,Connection}.DualError1(:,Stiff))./SimData{m,n,Connection}.SingleError1(:,Stiff); (SimData{m,n,Connection}.SingleError2(:,Stiff)-SimData{m,n,Connection}.DualError2(:,Stiff))./SimData{m,n,Connection}.SingleError2(:,Stiff)],'color',Color(Stiff,:),'marker','.','linestyle','none');

        [y_sim_fit,y_Sigma]=predict(SimErrorLME{m,n,Connection},dataset({[x,x.^2,zeros(length(x),1),log10(StiffArray(Stiff)*ones(length(x),1))],SimLMELabel{:}}));

        plot(x,y_sim_fit,'color',Color(Stiff,:),'linestyle','--','linewidth',1.5);
        fill([x;flipud(x)],[y_Sigma(:,1);flipud(y_Sigma(:,2))],Color(Stiff,:),'facealpha',0.2,'edgealpha',0);       

        % Data
        [y_data,y_Sigma]=predict(HHLME2,dataset({[x,x.^2,zeros(length(x),2),log10(StiffArray(Stiff)*ones(length(x),1))],LMELabel{:}}));
        plot(x,y_data,'color',Color(Stiff,:),'linewidth',1.5);
        fill([x;flipud(x)],[y_Sigma(:,1);flipud(y_Sigma(:,2))],Color(Stiff,:),'facealpha',0.2,'edgealpha',0);
        
        if Connection==2
            title('Corrected model');
        elseif Connection==3
            title('Naive model');
        elseif Connection==4
            title('Co-activity');
        elseif Connection==5
            title('Follow the better');
        end

        if Connection==5
            % Rescale axis for follow the better as it is out of view
            set(gca,'xtick',[-ScaleError,0,ScaleError/2],'ytick',-3:1:1);
            axis square; axis([-1,0.5,-3,1]*ScaleError); box on;               
        else
            set(gca,'xtick',[-ScaleError,0,ScaleError/2],'ytick',[-ScaleError/2,0,ScaleError/2]);
            axis image; axis([-1,0.5,-.5,0.5]*ScaleError); box on;
        end
        
        xlabel('Partner''s relative error (1-e_p/e)');
        ylabel('Improvement (1-e_c/e)');

        subplot(1,2,2); set(gca,'fontsize',SubplotFontSize); hold on;
        line([-ScaleError ScaleError],[0 0], 'Color', 'k', 'LineWidth', 1);
        line([0 0],[-.5 1]*ScaleEffort, 'Color', 'k', 'LineWidth', 1);
        
        % Simulation plot
        %plot([(SimData{m,n,Connection}.SingleError1(:,Stiff)-SimData{m,n,Connection}.SingleError2(:,Stiff))./SimData{m,n,Connection}.SingleError1(:,Stiff); (SimData{m,n,Connection}.SingleError2(:,Stiff)-SimData{m,n,Connection}.SingleError1(:,Stiff))./SimData{m,n,Connection}.SingleError2(:,Stiff)], [(SimData{m,n,Connection}.DualEffort1(:,Stiff)-SimData{m,n,Connection}.SingleEffort1(:,Stiff))./SimData{m,n,Connection}.SingleEffort1(:,Stiff); (SimData{m,n,Connection}.DualEffort2(:,Stiff)-SimData{m,n,Connection}.SingleEffort2(:,Stiff))./SimData{m,n,Connection}.SingleEffort2(:,Stiff)],'color',Color(Stiff,:),'marker','.','linestyle','none');

        [y_sim_fit,y_Sigma]=predict(SimEffortLME{m,n,Connection},dataset({[x,x.^2,x.^3,zeros(length(x),1),log10(StiffArray(Stiff)*ones(length(x),1))],SimEffortLMELabel{:}}));

        plot(x,y_sim_fit,'color',Color(Stiff,:),'linestyle','--','linewidth',1.5);
        fill([x;flipud(x)],[y_Sigma(:,1);flipud(y_Sigma(:,2))],Color(Stiff,:),'facealpha',0.2,'edgealpha',0);

        % Data
        [y_data,y_Sigma]=predict(HHEffortLME2,dataset({[x,x.^2,zeros(length(x),2),log10(StiffArray(Stiff)*ones(length(x),1))],LMELabel{:}}));
        plot(x,y_data,'color',Color(Stiff,:),'linewidth',1.5);
        fill([x;flipud(x)],[y_Sigma(:,1);flipud(y_Sigma(:,2))],Color(Stiff,:),'facealpha',0.2,'edgealpha',0);

        title(strcat('q=',num2str(SimData{1,1,1}.q(MinM)),',mu=',num2str(SimData{1,1,1}.mu(MinN)),',SSE=',num2str(MetricData(MinM,MinN,Connection))));
        xlabel('Partner''s relative error (1-e_p/e)');
        ylabel('Interaction effort (\alpha_c/\alpha-1)');

        set(gca,'xtick',[-ScaleError,0,ScaleError/2],'ytick',[-ScaleEffort/2,0,ScaleEffort]);
        axis square; axis([-ScaleError,ScaleError/2,-ScaleEffort/2,ScaleEffort]); box on;
    end
    
end



