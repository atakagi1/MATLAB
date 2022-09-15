function [Robot1,Robot2]= MechanismModel(Target,A,B,Q,L,SensoryNoise,SensoryDelay,Connection,k,d,dt,InternalModel)
    
    P1 = 0.001*diag(ones(1,size(A,1))); P2 = P1;
    Q1 = Q; Q2 = Q;
    L1 = L; L2 = L;

    InternalModel1 = InternalModel; InternalModel2 = InternalModel;
    IntModelDelay = SensoryDelay*ones(1,2);

    if Connection==0
        %%%%% Solo %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;
        H2 = H1;
        
        R1 = diag(SensoryNoise(1,1));
        R2 = diag(SensoryNoise(1,2));
        
        z1Noise = sqrt(SensoryNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(SensoryNoise(1,2))*randn(1,length(Target)-1);
    elseif Connection==1
        %%%%% Co-activity %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;           
        H2 = H1; 
        
        R1 = diag(SensoryNoise(1,1));
        R2 = diag(SensoryNoise(1,2));
        
        z1Noise = sqrt(SensoryNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(SensoryNoise(1,2))*randn(1,length(Target)-1);
    elseif Connection==2
        %%%%% Follow the better %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;          
        H2 = H1;        
        
        R1 = diag(SensoryNoise(1,1));
        R2 = diag(SensoryNoise(1,2));
        
        z1Noise = sqrt(SensoryNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(SensoryNoise(1,2))*randn(1,length(Target)-1);
    elseif Connection==3
        %%%% Multi-sensory integration %%%%
        H1 = zeros(2,size(A,1));
        H1(1,1) = 1;
        H1(2,1) = 1;
        H2 = H1;
        
        R1 = diag([SensoryNoise(1,1),SensoryNoise(1,2)]);
        R2 = diag([SensoryNoise(1,2),SensoryNoise(1,1)]);       
        
        z1Noise = [sqrt(SensoryNoise(1,1))*randn(1,length(Target)-1); sqrt(SensoryNoise(1,2))*randn(1,length(Target)-1)];
        z2Noise = [sqrt(SensoryNoise(1,2))*randn(1,length(Target)-1); sqrt(SensoryNoise(1,1))*randn(1,length(Target)-1)];
    elseif Connection==4
        %%%%% Interpersonal goal integration %%%%%
        H1 = zeros(2,size(A,1));
        H1(1,1) = 1;
        H1(2,1) = 1;
        H2 = H1;

        
        R1 = diag([SensoryNoise(1,1),SensoryNoise(1,2)]);
        R2 = diag([SensoryNoise(1,2),SensoryNoise(1,1)]);
        
        z1Noise = sqrt(R1)*randn(2,length(Target));
        z2Noise = sqrt(R2)*randn(2,length(Target));

        %%%%%% INTERNAL MODEL        
        A_1 = zeros(13);
        H_1 = zeros(7,size(A_1,1));
        H_1(1,1) = 1;
        H_1(2,2) = 1;
        H_1(3,3) = 1;
        H_1(4,7) = 1;
        H_1(5,8) = 1;
        H_1(6,9) = 1;
        H_1(7,10) = 1;
        
        H_2 = H_1;

        Q_1(1:3,1:3) = Q;
        Q_1(4:6,4:6) = Q;
        Q_1(7:9,7:9) = Q;
        Q_1(10,10) = 1;
        Q_1(11:13,11:13) = 10^-5*diag(ones(1,3));
        
        Q_2 = Q_1;
        
        P_1 = 0.001*diag(ones(1,size(A_1,1)));
        P_2 = P_1;
        R_1 = diag(SensoryNoise(1,1)*ones(1,size(H_1,1)));
        R_2 = diag(SensoryNoise(1,2)*ones(1,size(H_2,1)));
        xHat_1 = zeros(size(A_1,1),length(Target));
        xHat_2 = xHat_1;
        zNoise_1 = sqrt(R_1)*randn(size(H_1,1),length(Target)-1);
        zNoise_2 = sqrt(R_2)*randn(size(H_2,1),length(Target)-1);
        
        xHatProject_1 = zeros(size(A_1,1),length(Target));
        xHatProject_2 = xHatProject_1;
    
    elseif Connection==5
        %%%%% One-way (partner 2 is solo) %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;           
        H2 = H1; 
        
        R1 = diag(SensoryNoise(1,1));
        R2 = diag(SensoryNoise(1,2));
        
        z1Noise = sqrt(SensoryNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(SensoryNoise(1,2))*randn(1,length(Target)-1);
    end
    
    Robot1 = zeros(size(A,1),length(Target));
    Robot2 = zeros(size(A,1),length(Target));
   
    xHat1 = zeros(size(H1,2),length(Target));
    xHat2 = zeros(size(H2,2),length(Target));
    
    u1 = zeros(1,length(Target)); u2 = u1;
    F1 = zeros(1,length(Target)); F2 = F1;
    dF1 = zeros(1,length(Target)); dF2 = dF1;
    
    for i=2:length(Target)-2
        
        % Spring connection
        if Connection>0
            F1(1,i) = k*(Robot2(1,i)-Robot1(1,i))+d*(Robot2(2,i)-Robot1(2,i));
            F2(1,i) = -F1(1,i);
            
            dF1(1,i) = k*(Robot2(2,i)-Robot1(2,i))+d*(Robot2(3,i)-Robot1(3,i));
            dF2(1,i) = -dF1(1,i);
            
            if Connection==5
                F2(1,i) = 0;
                dF2(1,i) = 0;
            end
        end

        % Control policy + delay compensation
        if i>SensoryDelay
            xComp1 = xHat1(:,i-SensoryDelay);
            for s=1:SensoryDelay
                xComp1 = A*xComp1+B*(u1(1,i-SensoryDelay+s-1)+F1(1,i-SensoryDelay+s-1));
            end
            u1(1,i) = -L1*xComp1;    
        else
            u1(1,i) = 0;
        end
        if i>SensoryDelay
            xComp2 = xHat2(:,i-SensoryDelay);
            for s=1:SensoryDelay
                xComp2 = A*xComp2+B*(u2(1,i-SensoryDelay+s-1)+F2(1,i-SensoryDelay+s-1));
            end
            u2(1,i) = -L2*xComp2;
        else
            u2(1,i) = 0;
        end

        % Observation        
        if Connection==0
            if i>SensoryDelay
                z1 = Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay);
            else
                z1 = 0;
            end
            if i>SensoryDelay
                z2 = Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay);    
            else
                z2 = 0;
            end
        elseif Connection==4
            % Linearisation of nonlinear state
            A_1 = jaccsd(InternalModel1,xHat_1(:,i));
            A_2 = jaccsd(InternalModel2,xHat_2(:,i));
            
            if i>IntModelDelay(1,1)
                %%%%%%%%% INTERNAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                z_1 = [Robot1(1,i-IntModelDelay(1,1)); Robot1(2,i-IntModelDelay(1,1)); Robot1(3,i-IntModelDelay(1,1)); Target(1,i-IntModelDelay(1,1)); diff(Target(1,i-IntModelDelay(1,1):i+1-IntModelDelay(1,1)))/dt; diff(diff(Target(1,i-IntModelDelay(1,1):i+2-IntModelDelay(1,1)))/dt)/dt; F1(1,i-IntModelDelay(1,1))];

                %%% INTERNAL MODEL UPDATE %%%
                % Predict
                x_1=InternalModel1(xHat_1(:,i));
                P_1=A_1*P_1*A_1'+Q_1;
                K_1=P_1*H_1'/(H_1*P_1*H_1'+R_1);
                % Correct
                xHat_1(:,i+1)=x_1+K_1*(z_1+zNoise_1(:,i)-H_1*x_1);
                P_1=(eye(size(A_1))-K_1*H_1)*P_1;
                
                xHat_1(11:end,i+1) = max(xHat_1(11:end,i+1),0);
                
                
                % Project ahead to compensate for delay
                xIntTemp = xHat_1(:,i+1);
                for s=1:IntModelDelay(1,1)
                    % Predict
                    xIntTemp=InternalModel(xIntTemp);
                end
                xHatProject_1(:,i) = xIntTemp;
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % SENSORY MEASUREMENT
                z1 = [Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay); Robot1(1,i-SensoryDelay)-xHatProject_1(7,i)];
            else
                z1 = [0; 0];
            end
            if i>IntModelDelay(1,2)
                
                %%%%%%%%% INTERNAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                z_2 = [Robot2(1,i-IntModelDelay(1,2)); Robot2(2,i-IntModelDelay(1,2)); Robot2(3,i-IntModelDelay(1,2)); Target(1,i-IntModelDelay(1,2)); diff(Target(1,i-IntModelDelay(1,2):i+1-IntModelDelay(1,2)))/dt; diff(diff(Target(1,i-IntModelDelay(1,2):i+2-IntModelDelay(1,2)))/dt)/dt; F2(1,i-IntModelDelay(1,2))];

                %%% INTERNAL MODEL UPDATE %%%
                % Predict
                x_2=InternalModel2(xHat_2(:,i));
                P_2=A_2*P_2*A_2'+Q_2;
                K_2=P_2*H_2'/(H_2*P_2*H_2'+R_2);
                % Correct
                xHat_2(:,i+1)=x_2+K_2*(z_2+zNoise_2(:,i)-H_2*x_2);
                P_2=(eye(size(A_2))-K_2*H_2)*P_2;
                
                xHat_2(11:end,i+1) = max(xHat_2(11:end,i+1),0);
                
                % Project ahead to compensate for delay
                xIntTemp = xHat_2(:,i+1);
                for s=1:IntModelDelay(1,2)
                    % Predict
                    xIntTemp=InternalModel(xIntTemp);
                end
                xHatProject_2(:,i) = xIntTemp;
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % SENSORY MEASUREMENT
                z2 = [Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay); Robot2(1,i-SensoryDelay)-xHatProject_2(7,i)];
            else
                z2 = [0; 0];
            end
        elseif Connection==2
            if i>SensoryDelay
                if abs(Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay))<abs(Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay))
                    z1 = Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay);
                else
                    z1 = Robot1(1,i-SensoryDelay)-Robot2(1,i-SensoryDelay);
                end
                
                if abs(Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay))<abs(Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay))
                    z2 = Robot2(1,i-SensoryDelay)-Robot1(1,i-SensoryDelay);
                else
                    z2 = Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay);
                end
            else
                z1 = 0;
                z2 = 0;
            end
        elseif Connection==3
            if i>SensoryDelay
                z1 = [Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay); Robot1(1,i-SensoryDelay)-Robot2(1,i-SensoryDelay)];
            else
                z1 = [0;0];
            end
            if i>SensoryDelay
                z2 = [Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay); Robot2(1,i-SensoryDelay)-Robot1(1,i-SensoryDelay)];
            else
                z2 = [0;0];
            end
        elseif Connection==1||Connection==5
            % Co-activity and one-way
            if i>SensoryDelay
                z1 = Robot1(1,i-SensoryDelay)-Target(1,i-SensoryDelay);
            else
                z1 = 0;
            end
            if i>SensoryDelay
                z2 = Robot2(1,i-SensoryDelay)-Target(1,i-SensoryDelay);    
            else
                z2 = 0;
            end    
        end
        
        %%%% Partner 1 %%%%
        % Predict
        if i>SensoryDelay
            x1=A*xHat1(:,i-SensoryDelay)+B*(u1(1,i-SensoryDelay)+dF1(1,i-SensoryDelay));
            P1=A*P1*A'+Q1;
            K1=P1*H1'/(H1*P1*H1'+R1);
            z1 = z1+z1Noise(:,i);
            % Correct
            xHat1(:,i+1-SensoryDelay)=x1+K1*(z1-H1*x1);
            P1=(eye(length(K1*H1))-K1*H1)*P1;
        end

        %%%% Partner 2 %%%%
        % Predict
        if i>SensoryDelay
            x2=A*xHat2(:,i-SensoryDelay)+B*(u2(1,i-SensoryDelay)+dF2(1,i-SensoryDelay));
            P2=A*P2*A'+Q2;
            K2=P2*H2'/(H2*P2*H2'+R2);
            z2 = z2+z2Noise(:,i);
            % Correct
            xHat2(:,i+1-SensoryDelay)=x2+K2*(z2-H2*x2);
            P2=(eye(length(K2*H2))-K2*H2)*P2;
        end
        
        % Update actual system dynamics
        Robot1(:,i+1) = A*Robot1(:,i)+B*(u1(1,i)+dF1(1,i));
        Robot2(:,i+1) = A*Robot2(:,i)+B*(u2(1,i)+dF2(1,i));
    end
    
function A=jaccsd(fun,x)
% Jacobian through complex step differentiation
n=numel(x);
h=n*eps;
A = imag(fun(x(:,ones(1,n))+h*1i*eye(n)))/h;
 
    

