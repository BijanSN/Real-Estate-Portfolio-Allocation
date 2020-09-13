%If you have any trouble lunching this matlab file, please contact me, Bijan Sadat, at bijan.sadat@gmail.com 
%Thank you for reading.
%Requirements : Rmatlab2018 + financial toolbox functions
clear all, clc, close all

[ExchangeDate,SMI,SP500,CAC,FTSE,BONDUS,BONDCH,REUK,REFR,REASIA,REUS,REListedCH,REDirectCH] = importfile1('DATA_Total_Returns.xlsx','Feuil1',2,83);
Prices = [SMI,SP500,CAC,FTSE,BONDUS,BONDCH,REUK,REFR,REASIA,REUS,REListedCH,REDirectCH];
REUS= REUS(26:end)

%% Plots

%Stocks
plot(ExchangeDate,SP500)
hold on
plot(ExchangeDate,SMI)
plot(ExchangeDate,FTSE)
plot(ExchangeDate,CAC)

%Bonds
plot(ExchangeDate,BONDUS)
plot(ExchangeDate,BONDCH)

%Listed Real Estate
plot(ExchangeDate(26:end),REUS)
plot(ExchangeDate,REUK)
plot(ExchangeDate,REFR)
plot(ExchangeDate,REASIA)
plot(ExchangeDate,REListedCH)

%Direct Real estate
plot(ExchangeDate,REDirectCH)

legend('SP500','SMI','FTSE','CAC','BOND US','BOND CH','RE US','RE UK','RE FR','RE ASIA','Listed RE CH','Direct RE CH')
title('Nominal prices of Assets')
hold off

%% Plots ( returns)

Returns=(Prices(2:end,:)-Prices(1:end-1,:))./Prices(1:end-1,:)
Returns(1:26,10)=0;
Mean_Returns= sum(Returns,2,'omitnan')

Corr= corr(Returns)
% corrplot(Corr)
imagesc(Corr);
colormap(jet);
colorbar;

plot(ExchangeDate(2:end),Returns)% very volatile during 2008. 

TotalReturns=(Prices(end,:)-Prices(1,:))./Prices(1,:)
RE_RE_US= (REUS(end)-REUS(1))/REUS(1)
TotalReturns(10)= RE_RE_US

AnnualisedReturns=nthroot(TotalReturns,size(Returns,1)/4)-1 %81 (quaterly periods) /4 = number of yearly period : 20.25 year
%% Relative performence figures
% Let's use log returns( multiplicative)

Method='Continuous'
%Stocks
LR_SP500= tick2ret(SP500,'Method',Method)
LR_SMI= tick2ret(SMI,'Method',Method)
LR_FTSE=tick2ret(FTSE,'Method',Method)
LR_CAC=tick2ret(CAC,'Method',Method)

%Bonds
LR_BOND_US=  tick2ret(BONDUS,'Method',Method)
LR_BOND_CH=  tick2ret(BONDCH,'Method',Method)

%Listed Real Estate
LR_ListedRE_CH= tick2ret(REListedCH,'Method',Method)
LR_ListedRE_UK= tick2ret(REUK,'Method',Method)
LR_ListedRE_US= tick2ret(REUS,'Method',Method)
LR_ListedRE_ASIA= tick2ret(REASIA,'Method',Method)
LR_ListedRE_FR= tick2ret(REFR,'Method',Method)

%Direct Real Estate
LR_DirectRECH= tick2ret(REDirectCH,'Method',Method)

LReturns=[LR_BOND_US LR_BOND_CH  LR_ListedRE_CH LR_ListedRE_UK LR_ListedRE_FR LR_ListedRE_ASIA LR_SP500 LR_SMI LR_FTSE LR_CAC];
%LR_ListedRE_US
n=size(SP500,1); % Same timeframe
m=size(LReturns,2); % Number of assets
price =ones(n,m);
price(1,:) =1000;

for j=1:size(LReturns,2) %columns = each assets
    for i=1:(n-1)% row= each period
        if i==1
            price(i+1,j)= 1000+ price(i,j)*LReturns(i,j) %start at the same time to see performance
        else
            price(i+1,j)= price(i,j)+ price(i,j)*LReturns(i,j)
        end
    end
end

plot(ExchangeDate,price)
title('Relative performance of each assets' )

%% Box plot of Quaterly Returns

boxplot(Returns,'Labels',{'SMI','SP500','CAC','FTSE','BOND US','BOND CH','RE UK','RE FR','RE ASIA','RE US','RE Listed CH','RE Direct CH' } )

title('Boxplot of quaterly returns per assets')
xlabel('Assets')
ylabel('Returns')

%%

mret=mean(mean(Returns),'omitnan')% expected quaterly market returns
Amret= mean(mean(AnnualisedReturns,'omitnan')) %  annualised MARKET RETURNS

AssetCovar=cov(Returns) % covariance (unshrinked)
mrsk=mean(mean(AssetCovar)) % market risk

%rf
crsk=0;
cret=0.005; % rf = 0.5% annualised 
%% Setting up the portfolio p1 (need matlab 2018 + financial toolbox functions)

%p1 = Whole investment universe
AssetMean=mean(Returns)

p1=Portfolio;
p1=Portfolio(p1,'RiskFreeRate',cret)
p1=setAssetList(p1,'SMI','SP500','CAC','FTSE','BOND US','BOND CH','RE UK','RE FR','RE ASIA','RE US','RE Listed CH','RE Direct CH')
p1 = setAssetMoments(p1,AssetMean/sqrt(4),AssetCovar/sqrt(4));

p1 = setDefaultConstraints(p1);
pwgt = estimateFrontier(p1,100);

plotly=pwgt' %  > in plotly.

[prsk,pret] = estimatePortMoments(p1,pwgt);

%% Define specific portfolios among that universe p1 :

q = setBudget(p1, 0, 1);
qwgt = estimateFrontier(q,100);
[qrsk,qret] = estimatePortMoments(q,qwgt);

p1 = setInitPort(p1,1/p1.NumAssets);
[ersk,eret] = estimatePortMoments(p1,p1.InitPort); % e => EWP

portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'scatter', ersk, eret, {'EWP'}}, ...  
 	{'scatter', crsk, cret, {'Cash'}}, ...
    {'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'});

%% Creating portfolio 3 , without Real Estate : show benef of Real Estate
NoRE_AssetMean= AssetMean(1:6)/sqrt(4);
NoRE_Cov= cov(Returns(:,1:6))/sqrt(4);

p3=Portfolio;

p3=Portfolio(p3,'RiskFreeRate',cret)
p3=setAssetList(p3,'SMI','SP500','CAC','FTSE','BOND US','BOND CH');
p3 = setAssetMoments(p3,NoRE_AssetMean,NoRE_Cov);
p3 = setDefaultConstraints(p3);
pwgt3 = estimateFrontier(p3,100);
[prsk3,pret3] = estimatePortMoments(p3,pwgt3);

q3 = setBudget(p3, 0, 1);
qwgt3 = estimateFrontier(q3,100);
[qrsk3,qret3] = estimatePortMoments(q3,qwgt3);

p3 = setInitPort(p3,1/p3.NumAssets);

help estimateMaxSharpeRatio
[ersk3,eret3] = estimatePortMoments(p3,p3.InitPort) %EWP

portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk3, pret3}, ...
    {'line', qrsk3, qret3, [], [], 1}, ...
    {'scatter', ersk3, eret3, {'EWP Without Real Estate'}}, ...  
 	{'scatter', crsk, cret, {'Cash'}}, ... 
    {'scatter', sqrt(diag(p3.AssetCovar)), p3.AssetMean, p3.AssetList, '.r'});

% Suboptimal frontier
% EWP, for example, gets riskier without RE

%% Portfolio 12, everything but the DIRECT RE, to show benefits of having Direct RE !
NoDirectReturns= AssetMean(1:11)/sqrt(4);
NoDirectCov=cov(Returns(:,1:11))/sqrt(4);
p12=Portfolio;

p12=Portfolio(p12,'RiskFreeRate',cret)
p12=setAssetList(p12,'SMI','SP500','CAC','FTSE','BOND US','BOND CH','RE UK','RE FR','RE ASIA','RE US','RE Listed CH')
p12 = setAssetMoments(p12,NoDirectReturns,NoDirectCov);

p12 = setDefaultConstraints(p12);
pwgt12 = estimateFrontier(p12,100);

[prsk12,pret12] = estimatePortMoments(p12,pwgt12);

q12 = setBudget(p12, 0, 1);
qwgt12 = estimateFrontier(q12,100);
[qrsk12,qret12] = estimatePortMoments(q12,qwgt12);

portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk12, pret12}, ...
    {'line', qrsk12, qret12, [], [], 1}, ...
	{'scatter', crsk, cret, {'Cash'}}, ... 
    {'scatter', sqrt(diag(p12.AssetCovar)), p12.AssetMean, p12.AssetList, '.r'});

%% Plot WITH (p1)and WITHOUT RE (p3) together :
%Show inefficience if we do not include those assets
portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret, {'With Real Estate'},':b'}, ...
    {'line', qrsk, qret, [], [], 1}, ...  
    {'line', prsk3, pret3,{'Without Real Estate'}}, ...
    {'line', qrsk3, qret3, [], [], 1}, ...
    {'scatter', ersk, eret, {'EWP'}}, ...  
 	{'scatter', crsk, cret, {'Cash'}}, ...
    {'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'});

%% PLot all (p1) & all but direct (p12)

portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret,{'All assets'},':b'}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'line', prsk12, pret12,{'No Direct Real Estate'}}, ...
    {'line', qrsk12, qret12, [], [], 1}, ...
	{'scatter', crsk, cret, {'Cash'}}, ... 
    {'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'});

%% (!) low number of obs: correcton with another portfolio (schrinked)
cond(AssetCovar) %huge condition of the matrix : unstable

p2=Portfolio;
p2=Portfolio(p2,'RiskFreeRate',cret)
p2=setAssetList(p2,'S SMI','S SP500','S CAC','S FTSE','S BOND US','S BOND CH','S RE UK','S RE FR','S RE ASIA','S RE US','S RE Listed CH','S RE Direct CH')
[Mu_shr, Sigma_shr]=ShrinkLocationDispersion(Returns/sqrt(4)) % same ER, different Sigma !
p2 = setAssetMoments(p2,Mu_shr', Sigma_shr')

p2 = setDefaultConstraints(p2);
pwgt22 = estimateFrontier(p2,100);

[prsk22,pret22] = estimatePortMoments(p2,pwgt22);
plotFrontier(p2) 

q2 = setBudget(p2, 0, 1);
qwgt22 = estimateFrontier(q2,100);
[qrsk22,qret22] = estimatePortMoments(q2,qwgt22);
% plotFrontier(q)

p2= setInitPort(p2,1/p2.NumAssets);
[ersk22,eret22] = estimatePortMoments(p2,p2.InitPort); % e => EWP


portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret,{'Unschrinked'},':b'}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'line', prsk22, pret22,{'Schrinked'}}, ...
    {'line', qrsk22, qret22, [], [], 1}, ...
    {'scatter', ersk, eret, {'EWP UN S'}}, ... 
 	{'scatter', ersk22, eret22, {'EWP S'}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ...
    {'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'},...
    {'scatter', sqrt(diag(p2.AssetCovar)), p2.AssetMean, p2.AssetList, '.r'});


%% Using OPP2 constraints on P1
% max 50% Equity, max 30% RE

%4 first =  Stocks
% then 2 Bonds
% then 6 Real Estate
Group = zeros(1,12);
Group(1,1:4)=1; %Stocks
Group(2,7:12)=1; %Real estate

GroupMin = [0 0];
GroupMax = [0.5 0.3];

[A,b] = pcglims(Group,GroupMin,GroupMax)

NumPorts=100;
p6=Portfolio;

p6=Portfolio(p6,'RiskFreeRate',cret)
p6=setAssetList(p6,'SMI','SP500','CAC','FTSE','BOND US','BOND CH','RE UK','RE FR','RE ASIA','RE US','RE Listed CH','RE Direct CH')
p6 = setAssetMoments(p6,AssetMean/sqrt(4),AssetCovar/sqrt(4));

p6 = setBounds(p6, 0, 1);
p6 = setBudget(p6, 1, 1);

p6.AInequality = A
p6.bInequality =b

PortWts6 = estimateFrontier(p6, NumPorts);

Weigths_constrainted= PortWts6'
[PortRisk6, PortReturn6] = estimatePortMoments(p6, PortWts6);

plotly6=PortWts6' %plotly

[prsk6,pret6] = estimatePortMoments(p6,PortWts6);
q = setBudget(p6, 0, 1);
qwgt6 = estimateFrontier(q,100);
[qrsk6,qret6] = estimatePortMoments(q,PortWts6);
p6 = setInitPort(p6,1/p6.NumAssets);
[ersk6,eret6] = estimatePortMoments(p6,p6.InitPort); % e => EWP

%frontière avec contraintes et sans, together :
close all
%pf1
portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret,{'Unconstrainted'},':b'}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'line', prsk6, pret6,{'OPP2 Constrained'}}, ...
    {'line', qrsk6, qret6, [], [], 1}, ... %     {'scatter', arsk, aret, {sprintf('%g%% Return',100*TargetReturn)}}, ...%     {'scatter', brsk, bret, {sprintf('%g%% Risk',100*TargetRisk)}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ... %{'scatter', Srsk, Sret, {'Constraint no Real Estate'}}, ...
    {'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'});

%% Sans RE FR,car outlier : 11 assets now.

ReturnsNOFR= [Returns(:,1:7) Returns(:,9:end)];
AssetMeannoFR= mean(ReturnsNOFR);
AssetCovarNOFR=cov(ReturnsNOFR) % covariance 

p=Portfolio;

p=Portfolio(p,'RiskFreeRate',cret)
p=setAssetList(p,'SMI','SP500','CAC','FTSE','BOND US','BOND CH','RE UK','RE ASIA','RE US','RE Listed CH','RE Direct CH')
p = setAssetMoments(p,AssetMeannoFR/sqrt(4),AssetCovarNOFR/sqrt(4));

p = setBounds(p, 0, 1);
p = setBudget(p, 1, 1);

Group = zeros(2,11);
Group(1,1:4)=1; %Stocks 
Group(2,7:11)=1; %Real estate

GroupMin = [0 0];
GroupMax = [0.5 0.3];

[A,b] = pcglims(Group,GroupMin,GroupMax)

p.AInequality = A
p.bInequality =b

PortWtsp = estimateFrontier(p, NumPorts);

Weigths_constrainted_noFR= PortWtsp'
[PortRisk, PortReturn] = estimatePortMoments(p, PortWtsp)

%% same, 11 assets, no constraints

p=Portfolio; %reset : no constraints
p=setAssetList(p,'SMI','SP500','CAC','FTSE','BOND US','BOND CH','RE UK','RE ASIA','RE US','RE Listed CH','RE Direct CH')

p=Portfolio(p,'RiskFreeRate',cret)
p = setAssetMoments(p,AssetMeannoFR/sqrt(4),AssetCovarNOFR/sqrt(4));

p = setBounds(p, 0, 1);
p = setBudget(p, 1, 1);

PortWtsp123 = estimateFrontier(p, NumPorts);

Weigths_unconstrainted_noFR= PortWtsp123'
[PortRisk, PortReturn] = estimatePortMoments(p, PortWtsp)

%% Value at risk of EWP (p1= all asset unschrinked) 

p1 = setInitPort(p1,1/p1.NumAssets);
z_95=norminv(.95,0,1); %Inverse normal cumulative distribution function giving the 95% quantile
z_99=norminv(.99,0,1);

a=p1.InitPort
mu=AssetMean;
SIGMA=AssetCovar;
z=-AssetMean*a; %This is the negative of the returns(=losses) weighted by their portfolio weight (a)
alpha=0.1;

GaussianVaR95=-(z-(a'*SIGMA*a)^(1/2)*z_95)
GaussianVaR=-(z-(a'*SIGMA*a)^(1/2)*z_99) % Value-at-Risk of all assets 

opt=optimset('Display','iter'); %Optimization parameters
guess=GaussianVaR; %guess for the numerical minimization procedure.
NonParamVaR=fminsearch(@obj_fct, guess, opt, z) % same as the gaussian one

%% Value at risk of EWP (p3= no Real Estate) 

p3 = setInitPort(p3,1/p3.NumAssets);
z_95=norminv(.95,0,1); %Inverse normal cumulative distribution function giving the 95% quantile
z_99=norminv(.99,0,1);

a=p3.InitPort

mu=NoRE_AssetMean;
SIGMA=NoRE_Cov;
z=-mu*a; %This is the negative of the returns(=losses) weighted by their portfolio weight (a)
alpha=0.1;

GaussianVaR95_NORE=-(z-(a'*SIGMA*a)^(1/2)*z_95)
GaussianVaR=-(z-(a'*SIGMA*a)^(1/2)*z_99) %Lower value at risk, because Real estate = riskier than stocks here

%% final table

plotFrontier(p1)
w1=estimateMaxSharpeRatio(p1)
[risk, ret] = estimatePortMoments(p1, w1);
hold on
plot(risk, ret,'or') % all assets
w2=estimateMaxSharpeRatio(p12)
[risk2, ret2] = estimatePortMoments(p12, w2);
plot(risk2, ret2,'ob') % no DIRECT REAL ESTATE
w3=estimateMaxSharpeRatio(p2)
[risk3, ret3] = estimatePortMoments(p2, w3);
plot(risk3, ret3,'og') % schrinked
w4=estimateMaxSharpeRatio(p3) 
[risk4, ret4] = estimatePortMoments(p3, w4);
plot(risk4, ret4,'ok')% No Real Estate
w5=estimateMaxSharpeRatio(p6)
[risk5, ret5] = estimatePortMoments(p6, w5);
plot(risk5, ret5,'o') % OPP2 Constrained
legend('EWP all Assets', 'Efficient Frontier','All Assets', 'No direct RE', 'Schrinked', 'No Real Estate', 'OPP2 Constrained')