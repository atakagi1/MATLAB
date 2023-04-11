function [Robot1,Robot2,uRecord1,uRecord2]= compliance_model2(Target,A,B,Q,L,TaskNoise,Connection,k,ErrorToNoiseParam,NoiseToErrorParam)
       
    P1 = diag(ones(1,size(A,1)))*0.001;
    P2 = P1;

    if Connection==0
        %%%%% Solo %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;
            
        H2 = H1;
        
        R1 = diag(TaskNoise(1,1));
        R2 = diag(TaskNoise(1,2));
        
        z1Noise = sqrt(TaskNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(TaskNoise(1,2))*randn(1,length(Target)-1);
    elseif Connection==1
        %%%%% RELAXATION %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;
            
        H2 = H1;
        
        R1 = diag(TaskNoise(1,1));
        R2 = diag(TaskNoise(1,2));
        
        z1Noise = sqrt(TaskNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(TaskNoise(1,2))*randn(1,length(Target)-1);
    elseif Connection==2||Connection==3
        %%%%% Two way %%%%%
        H1 = zeros(2,size(A,1));
        H1(1,1) = 1;
        H1(2,1) = 1;

        H2 = H1;
        
        % Set connection strength and additional error from compliance
        % (this comes from a separate experiment)
        switch k
            case 0.3*180/pi
                ErrorAdd = 0.4372;
            case 0.03*180/pi
                ErrorAdd = 3.0268;
            case 0.005*180/pi
                ErrorAdd = 9.7755;
        end
        
        if Connection==3
            ErrorAdd = 0;
        end

        R1 = diag([TaskNoise(1,1),polyval(ErrorToNoiseParam,ErrorAdd+polyval(NoiseToErrorParam,sqrt(TaskNoise(1,2)*(180/pi)^2)))^2*(pi/180)^2]);
        R2 = diag([TaskNoise(1,2),polyval(ErrorToNoiseParam,ErrorAdd+polyval(NoiseToErrorParam,sqrt(TaskNoise(1,1)*(180/pi)^2)))^2*(pi/180)^2]);
        
        z1Noise = sqrt(R1)*randn(size(R1,1),length(Target)-1);
        z2Noise = sqrt(R2)*randn(size(R2,1),length(Target)-1);
    elseif Connection==4
        %%%%% Co-activity %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;
            
        H2 = H1;
        
        R1 = diag(TaskNoise(1,1));
        R2 = diag(TaskNoise(1,2));
        
        z1Noise = sqrt(TaskNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(TaskNoise(1,2))*randn(1,length(Target)-1);
    elseif Connection==5
        %%%%% Follow the better %%%%%
        H1 = zeros(1,size(A,1));
        H1(1,1) = 1;
            
        H2 = H1;
        
        R1 = diag(TaskNoise(1,1));
        R2 = diag(TaskNoise(1,2));
        
        z1Noise = sqrt(TaskNoise(1,1))*randn(1,length(Target)-1);
        z2Noise = sqrt(TaskNoise(1,2))*randn(1,length(Target)-1);
    end
    
    Robot1 = zeros(size(A,1),length(Target));
    Robot2 = zeros(size(A,1),length(Target));
   
    xHat1 = Robot1;
    xHat2 = Robot2;
    
    F1 = 0; F2 = 0;
    
    uRecord1 = zeros(1,length(Target));
    uRecord2 = zeros(1,length(Target));
    
    for i=1:length(Target)-1
        
        % Spring connection
        if Connection~=0
            HapticForce = k*(Robot2(1,i)-Robot1(1,i))+(0.005*180/pi)/10*(Robot2(2,i)-Robot1(2,i));
            
            F1 = HapticForce;
            F2 = -HapticForce;
        end

        % Control policy
        if i>1
            u1 = -L*xHat1(:,i);
            u2 = -L*xHat2(:,i);
        else
            u1 = 0;
            u2 = 0;
        end

        % "Real" system
        Robot1(:,i+1) = A*Robot1(:,i)+B*(u1+F1);
        Robot2(:,i+1) = A*Robot2(:,i)+B*(u2+F2);
        
        % Observation        
        if Connection==0            
            z1 = Robot1(1,i+1)-Target(1,i+1);
            z2 = Robot2(1,i+1)-Target(1,i+1);
        elseif Connection==1
            z1 = Robot1(1,i+1)-Target(1,i+1);
            z2 = Robot2(1,i+1)-Target(1,i+1);
            
            e1 = abs(Robot1(1,i+1)-Target(1,i+1));
            e2 = abs(Robot2(1,i+1)-Target(1,i+1));

            if e1>=e2
                u1 = 0;
            else
                u2 = 0;
            end
        elseif Connection==2||Connection==3
            z1 = [Robot1(1,i+1)-Target(1,i+1); Robot1(1,i+1)-Target(1,i+1)];
            z2 = [Robot2(1,i+1)-Target(1,i+1); Robot2(1,i+1)-Target(1,i+1)];
        elseif Connection==4
            z1 = Robot1(1,i+1)-Target(1,i+1);
            z2 = Robot2(1,i+1)-Target(1,i+1);
        elseif Connection==5
            e1 = abs(Robot1(1,i+1)-Target(1,i+1));
            e2 = abs(Robot2(1,i+1)-Target(1,i+1));
            if e1>=e2
                z1 = Robot1(1,i+1)-Robot2(1,i+1);
                z2 = Robot2(1,i+1)-Target(1,i+1);
            else
                z1 = Robot1(1,i+1)-Target(1,i+1);
                z2 = Robot2(1,i+1)-Robot1(1,i+1);
            end
        end  
        
        %%%% Partner 1 %%%%
        % Predict
        x1=A*xHat1(:,i)+B*(u1+F1);
        P1=A*P1*A'+Q;
        K1=P1*H1'/(H1*P1*H1'+R1);
        z1 = z1+z1Noise(:,i);
        
        % Correct
        xHat1(:,i+1)=x1+K1*(z1-H1*x1);
        P1=(eye(length(K1*H1))-K1*H1)*P1;    
        
        %%%% Partner 2 %%%%
        % Predict
        x2=A*xHat2(:,i)+B*(u2+F2);
        P2=A*P2*A'+Q;
        K2=P2*H2'/(H2*P2*H2'+R2);
        z2 = z2+z2Noise(:,i);
        
        % Correct
        xHat2(:,i+1)=x2+K2*(z2-H2*x2);
        P2=(eye(length(K2*H2))-K2*H2)*P2;
        

        uRecord1(1,i) = u1;
        uRecord2(1,i) = u2;
    end
    
    % Convert system from radians to degrees
    Robot1 = Robot1*180/pi;
    Robot2 = Robot2*180/pi;
    
    
    
    

