function [ConfBand,beta,ModelFunction,h,hFill] = HACPlot(x,y,PlotColor,Order,opts)

if nargin<=4
    LineStyle = '-';
    LineWidth = 2;
else
    if opts==0
        LineStyle = '-';
    else
        LineStyle = '--';
    end
    LineWidth = 3;
end
XVector = transpose(min(x):0.01:max(x));

% Model function
ModelFunction = @(b,x)(polyval(fliplr(b),x));

% Non-linear regression
Mdl = fitlm(x,y,strcat('poly',num2str(Order)));
opts = statset('nlinfit');
opts.RobustWgtFun = 'huber';
[beta,nlresid,J] = nlinfit(x,y,ModelFunction,ones(1,Order+1),opts);

%NWEstParamCov = hac(Mdl);
%[~,NWfitcb] = nlpredci(ModelFunction,XVector,beta,nlresid,'Covar',NWEstParamCov); % Corrected margin of error

[~,NWfitcb] = nlpredci(ModelFunction,XVector,beta,nlresid,'jacobian',J,'SimOpt','off'); % Corrected margin of error
ConfBand = [ModelFunction(beta,XVector)-NWfitcb, ModelFunction(beta,XVector)+NWfitcb]; % Corrected confidence bands


h=plot(XVector,ModelFunction(beta,XVector),'color',PlotColor,'linewidth',LineWidth,'linestyle',LineStyle);

if nargin<=4
    hFill=fill([XVector;flipud(XVector)],[ConfBand(:,1);flipud(ConfBand(:,2))],PlotColor);
    set(hFill,'facealpha',0.2,'edgecolor','none','edgealpha',0);
end
