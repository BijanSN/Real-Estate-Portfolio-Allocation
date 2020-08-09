clear all, clc

%% loading data

%USA
[ExchangeDate,SP500] = importfile('sp500_monthly .xlsx','Sheet 1',15,114);
[ExchangeDate3,RE_US] = importfile4('re_us.xlsx','Sheet 1',15,114);
[Date1,BOND_US] = importfile8('bond us government .xlsx','Worksheet',7,106);

%CH
[ExchangeDate2,SMI] = importfile1('SMI_monthly.xlsx','Sheet 1',26,125);
[ExchangeDate1,ListedRE_CH] = importfile2('RE_CH_TRUE.xlsx','Sheet 1',15,114);
[Date,BOND_CH] = importfile7('bond CH_gov.xlsx','Worksheet',8,107);

% [AppartementsenPPE1,Villas1] = importfile9('Direct RE Switzerland.xlsx','Feuille 1',8,106)

RECHTRUE = importfile10('RE_CH_TRUE.xlsx','Sheet 1',15,114); %  à verif si bon file.

%UK
[ExchangeDate7,FTSE] = importfile3('FTSE_monthly .xlsx','Sheet 1',24,123);
[ExchangeDate8,RE_UK] = importfile5('re_uk.xlsx','Sheet 1',15,114);

%Date==ExchangeDate2


%semillogy possible
plot(ExchangeDate,BOND_US)
hold on
plot(ExchangeDate,BOND_CH)
plot(ExchangeDate,RE_US)
plot(ExchangeDate,ListedRE_CH)
plot(ExchangeDate,RE_UK)
plot(ExchangeDate, SP500)
plot(ExchangeDate,SMI)
plot(ExchangeDate,FTSE)
% plot(ExchangeDate,AppartementsenPPE1)
% plot(ExchangeDate,Villas1)
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


% plot(R_BOND_CH,'k-')
% hold on
% plot(R_ListedRE_CH, 'g-')
% plot(R_SMI,'r-')
% 
% std(R_BOND_CH)
% std(R_ListedRE_CH)
% std(R_SMI)
% hold off

Returns=[R_BOND_US R_BOND_CH R_ListedRE_US R_ListedRE_CH R_ListedRE_UK R_SP500 R_SMI R_FTSE];
plot(Returns)
title('All returns')

% to annualise /?\










 % /!\assumption I made ( to correct ! ) :  our "market return" =  Average returns of all assets (1/N)

AssetMean=mean(Returns)
mret=mean(AssetMean); %market returns (as the mean of every asset over the whole timeframe)
AssetCovar=cov(Returns) % covariance (unshrinked)
mrsk=mean(mean(AssetCovar)) % market risk


%% using fmincon 
help optimoptions
o=optimoptions(fmincon, 'Algorithm','interior-point');
o.TolFun =1e-8;






%% rf
crsk=0;
cret=-0.004; % TO CHANGE WRT DATA

%% Boxplot + Summary statistics

boxplot(Returns,'Labels',{'BOND US','BOND CH','ListedRE US','ListedRE CH','ListedRE UK','SP500','SMI','FTSE'} )
%'Notch','on'
title('Boxplot of Returns per assets')
xlabel('Assets')
ylabel('Returns')

% Faire un tableau limite avec tout les assets

% E_BOND_US= mean(R_BOND_US)
% S_BOND_US= std(R_BOND_US)
% plot(S_BOND_US, S_BOND_US/E_BOND_US)
% xlabel('Assets')
% ylabel('Returns')

% std(R_BOND_CH)
% std(R_ListedRE_US)
% std(R_ListedRE_CH)
% std(R_ListedRE_UK)
% std(R_SP500)
% std(R_SMI)
% std(R_FTSE)

%% Setting up the portfolio

p1=Portfolio;
p1=setAssetList(p1,'BOND US','BOND CH','ListedRE US','ListedRE CH','ListedRE UK','SP500','SMI','FTSE')
p1 = setAssetMoments(p1,AssetMean,AssetCovar);

% Asset mean as historical mean of returns
% Covar = basic covar

%% 

p1 = setDefaultConstraints(p1);
pwgt = estimateFrontier(p1,100);

plotly=pwgt' % les weigths transposées=> dans plotly.
% https://plot.ly/create/?fid=BijanSN:3&fid=BijanSN:2

[prsk,pret] = estimatePortMoments(p1,pwgt);
% plotFrontier(p)

q = setBudget(p1, 0, 1);
qwgt = estimateFrontier(q,100);
[qrsk,qret] = estimatePortMoments(q,qwgt);


p1 = setInitPort(p1,1/p1.NumAssets);
[ersk,eret] = estimatePortMoments(p1,p1.InitPort); %EWP

% si contraintes de pas investir dans asset 5 par ex: 

% p.InitPort = [0.25 0.25 0.25 0.25 0]';
% [Srsk,Sret] = estimatePortMoments(p,p.InitPort); %PF sans asset 5

% TargetReturn = 1.3*eret;            % input target annualized return and risk here
% TargetRisk = 0.8*ersk;
% 
% awgt = estimateFrontierByReturn(p,TargetReturn/12);
% [arsk,aret] = estimatePortMoments(p,awgt);
% 
% bwgt = estimateFrontierByRisk(p,TargetRisk/sqrt(12));
% [brsk,bret] = estimatePortMoments(p,bwgt);


clf;
subplot(1,2,1)
portfolioexamples_plot('Asset Risks and Returns', ... 
    {'line', prsk, pret}, ...
    {'line', qrsk, qret, [], [], 1}, ...
    {'scatter', ersk, eret, {'EWP'}}, ...%   {'scatter', arsk, aret, {sprintf('%g%% Return',100*TargetReturn)}}, ...% 	{'scatter', brsk, bret, {sprintf('%g%% Risk',100*TargetRisk)}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ... % {'scatter', Srsk, Sret, {'Constraint no asset 5'}}, ...
	{'scatter', sqrt(diag(p1.AssetCovar)), p1.AssetMean, p1.AssetList, '.r'});


%% Try number 2, with another portfolio (schrinked)

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
    {'scatter', ersk, eret, {'EWP'}}, ...%   {'scatter', arsk, aret, {sprintf('%g%% Return',100*TargetReturn)}}, ...% 	{'scatter', brsk, bret, {sprintf('%g%% Risk',100*TargetRisk)}}, ...
 	{'scatter', crsk, cret, {'Cash'}}, ... % {'scatter', Srsk, Sret, {'Constraint no asset 5'}}, ...
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
%0.008

% GaussianES=mu+SIGMA*normpdf(z_95,0,1)/alpha
%% With constraints: no real estate (listed)

p1.InitPort = [0.2 0.2 0 0 0 0.2 0.2 0.2]';% no real estate (listed)
a=p1.InitPort

GaussianVaR=-(z-(a'*SIGMA*a)^(1/2)*z_95)

opt=optimset('Display','iter'); %Optimization parameters
guess=GaussianVaR; %guess for the numerical minimization procedure.
NonParamVaR=fminsearch(@obj_fct, guess, opt, z)

%lower var than unconstrained ! => not good.

% GaussianES=mu+SIGMA*normpdf(z_95,0,1)/alpha

