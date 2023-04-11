function [yFit,NWfitcb,ModelFunction,h,hFill] = HACPlot(x,y,PlotColor,Order,LineWidth,FaceAlpha,LineStyle)

XVector = transpose(-4:0.1:4);%transpose(min(x):0.1:max(x));

% Model function
ModelFunction = @(b,x)(polyval(fliplr(b),x));

% Non-linear regression
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[beta,nlresid,J] = nlinfit(x,y,ModelFunction,ones(1,Order+1),opts);

yFit = ModelFunction(beta,XVector);
[~,NWfitcb] = nlpredci(ModelFunction,XVector,beta,nlresid,'Jacobian',J); % Corrected margin of error
ConfBand = [yFit-NWfitcb, yFit+NWfitcb]; % Corrected confidence bands


h=plot(XVector,yFit,'color',PlotColor,'linewidth',LineWidth,'linestyle',LineStyle);
hFill=fill([XVector;flipud(XVector)],[ConfBand(:,1);flipud(ConfBand(:,2))],PlotColor);
set(hFill,'facealpha',0.2,'edgecolor','none','edgealpha',0);
