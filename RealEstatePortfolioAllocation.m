clear all, clc

%% loading data

%USA
[ExchangeDate,SP500] = importfile('sp500_monthly.xlsx','Sheet 1',15,114);
[ExchangeDate3,RE_US] = importfile4('re_us.xlsx','Sheet 1',15,114);
[Date1,BOND_US] = importfile8('bond us government.xlsx','Worksheet',7,106);

%CH
[ExchangeDate2,SMI] = importfile1('SMI_monthly.xlsx','Sheet 1',26,125);
[ExchangeDate1,ListedRE_CH] = importfile2('RE_CH_TRUE.xlsx','Sheet 1',15,114);
[Date,BOND_CH] = importfile7('bond CH_gov.xlsx','Worksheet',8,107);

%UK
[ExchangeDate7,FTSE] = importfile3('FTSE_monthly.xlsx','Sheet 1',24,123);
[ExchangeDate8,RE_UK] = importfile5('re_uk.xlsx','Sheet 1',15,114);

%%  
plot(ExchangeDate,BOND_US)
hold on
plot(ExchangeDate,BOND_CH)
plot(ExchangeDate,RE_US)
plot(ExchangeDate,ListedRE_CH)
plot(ExchangeDate,RE_UK)
plot(ExchangeDate, SP500)
plot(ExchangeDate,SMI)
plot(ExchangeDate,FTSE)

legend('BOND US','BOND CH','ListedRE US','ListedRE CH','ListedRE UK','SP500','SMI','FTSE')
hold off
%% Computation of mean returns & covar

%USA

R_SP500= tick2ret(SP500)
R_ListedRE_US= tick2ret(RE_US)
R_BOND_US=  tick2ret(BOND_US)
%CH
R_SMI= tick2ret(SMI)
R_ListedRE_CH= tick2ret(ListedRE_CH)
R_BOND_CH=  tick2ret(BOND_CH)

%UK
R_FTSE=tick2ret(FTSE)
R_ListedRE_UK= tick2ret(RE_UK)

Returns=[R_BOND_US R_BOND_CH R_ListedRE_US R_ListedRE_CH R_ListedRE_UK R_SP500 R_SMI R_FTSE];
plot(Returns)
title('All returns')

AssetMean=mean(Returns)
mret=mean(AssetMean); %market returns (as the mean of every asset over the whole timeframe)
AssetCovar=cov(Returns) % covariance (unshrinked)
mrsk=mean(mean(AssetCovar)) % market risk

%o=optimoptions(fmincon, 'Algorithm','interior-point');
o.TolFun =1e-8;

crsk=0;
cret=-0.004;
%% Relative performence figure
% let's use log returns ( multiplicative)

Method='Continuous'

LR_SP500= tick2ret(flipud(SP500),'Method',Method)
LR_ListedRE_US= tick2ret(flipud(RE_US),'Method',Method)
LR_BOND_US=  tick2ret(flipud(BOND_US),'Method',Method)

%CH
LR_SMI= tick2ret(flipud(SMI),'Method',Method)
LR_ListedRE_CH= tick2ret(flipud(ListedRE_CH),'Method',Method)
LR_BOND_CH=  tick2ret(flipud(BOND_CH),'Method',Method)

%UK
LR_FTSE=tick2ret(flipud(FTSE),'Method',Method)
LR_ListedRE_UK= tick2ret(flipud(RE_UK),'Method',Method)

LReturns=[LR_BOND_US LR_BOND_CH LR_ListedRE_US LR_ListedRE_CH LR_ListedRE_UK LR_SP500 LR_SMI LR_FTSE];

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

plot(flipud(ExchangeDate),price)
title('Relative performance of each assets' )

%% Boxplot + Summary statistics

boxplot(Returns,'Labels',{'BOND US','BOND CH','ListedRE US','ListedRE CH','ListedRE UK','SP500','SMI','FTSE'} )
title('Boxplot of Returns per assets')
xlabel('Assets')
ylabel('Returns')

%% Setting up the portfolio

p1=Portfolio;
p1=setAssetList(p1,'BOND US','BOND CH','ListedRE US','ListedRE CH','ListedRE UK','SP500','SMI','FTSE')
p1 = setAssetMoments(p1,AssetMean,AssetCovar);

%% 

p1 = setDefaultConstraints(p1);
pwgt = estimateFrontier(p1,100);

plotly=pwgt'
% https://plot.ly/create/?fid=BijanSN:3&fid=BijanSN:2

[prsk,pret] = estimatePortMoments(p1,pwgt);

q = setBudget(p1, 0, 1);
qwgt = estimateFrontier(q,100);
[qrsk,qret] = estimatePortMoments(q,qwgt);

p1 = setInitPort(p1,1/p1.NumAssets);
[ersk,eret] = estimatePortMoments(p1,p1.InitPort); %EWP

clf;
subplot(1,2,1)
portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'scatter', ersk, eret, {'EWP'}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ... 
	{'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'});


%% Portfolio without Real Estate :

p3=Portfolio;
p3=setAssetList(p3,'BOND US','BOND CH','SP500','SMI','FTSE')

AssetMean3= AssetMean
AssetMean3(:,3:5)=[]

AssetCovar3=AssetCovar
AssetCovar3(:,3:5)=[]
AssetCovar3(3:5,:)=[]

p3 = setAssetMoments(p3,AssetMean3,AssetCovar3);

p3 = setDefaultConstraints(p3);
pwgt3 = estimateFrontier(p3,100);

plotly=pwgt3'

[prsk3,pret3] = estimatePortMoments(p3,pwgt3);

q3 = setBudget(p3, 0, 1);
qwgt3 = estimateFrontier(q3,100);
[qrsk3,qret3] = estimatePortMoments(q3,qwgt3);

p3 = setInitPort(p3,1/p3.NumAssets);
[ersk3,eret3] = estimatePortMoments(p3,p3.InitPort); %EWP

subplot(1,2,2)
portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk3, pret3}, ...
    {'line', qrsk3, qret3, [], [], 1}, ...
    {'scatter', ersk3, eret3, {'EWP'}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ... 
	{'scatter', sqrt(diag(p3.AssetCovar)), p3.AssetMean, p3.AssetList, '.r'});


%% Portfolio number 2,Covar matrix schrinked

p2=Portfolio;
p2=setAssetList(p2,'BOND US','BOND CH','ListedRE US','ListedRE CH','ListedRE UK','SP500','SMI','FTSE')


[Mu_shr, Sigma_shr]=ShrinkLocationDispersion(Returns) % same ER, different Sigma !
%before :
%p = setAssetMoments(p,AssetMean,AssetCovar);
p2 = setAssetMoments(p2,Mu_shr', Sigma_shr')


subplot(1,2,2)
portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'scatter', ersk, eret, {'EWP'}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ...
	{'scatter', sqrt(diag(p2.AssetCovar)), p2.AssetMean, p2.AssetList, '.r'});

%% various stats on EWP ( schrinked) 

p1 = setInitPort(p1,1/p1.NumAssets);
z_95=norminv(.95,0,1); %Inverse normal cumulative distribution function giving the 95% quantile
z_99=norminv(.99,0,1);

a=p1.InitPort
mu=AssetMean;
SIGMA=AssetCovar;
z=-AssetMean*a; %This is the negative of the returns(=losses) weighted by their portfolio weight (a)
alpha=0.1;

GaussianVaR=-(z-(a'*SIGMA*a)^(1/2)*z_95)

opt=optimset('Display','iter'); %Optimization parameters
guess=GaussianVaR; %guess for the numerical minimization procedure.
NonParamVaR=fminsearch(@obj_fct, guess, opt, z)

%% With constraints: no real estate (listed)

p1.InitPort = [0.2 0.2 0 0 0 0.2 0.2 0.2]';% no real estate (listed)
a=p1.InitPort

GaussianVaR=-(z-(a'*SIGMA*a)^(1/2)*z_95)

opt=optimset('Display','iter'); %Optimization parameters
guess=GaussianVaR; %guess for the numerical minimization procedure.
NonParamVaR=fminsearch(@obj_fct, guess, opt, z)
