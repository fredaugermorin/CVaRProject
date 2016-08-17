%% Replicate Credit Manager Credit VaR simulation Procedure provided in their excel spreadsheet
clear all;clc;
%% Load preprocessed data
load('data.met.mat')

%% User inputs
decompositionMethod = 'cholesky'; %Choose cholesky or SVD
nbSim = 1000000; %number of simulations
meanRR = 0.58; %mean of default recovery rates
stdRR = 0.23; %STANDARD DEVIATION of default recovery rates
recoveryClaim = 1014426.23; %Recovery claim

%% Perform several computations to obtain weights on synthetic factors.
% First obtain implied idiosyncratic weights
load('data.met.mat')

idioW = zeros(length(R2),length(R2));
stdDev = zeros(size(R2));
for ii = 1:length(R2)
    idioW(ii,ii) = sqrt((1/R2(ii)-1)*relativeW(ii,:)*varCov*relativeW(ii,:)');
    stdDev(ii) = sqrt(relativeW(ii,:)*varCov*relativeW(ii,:)' +  idioW(ii,ii)^2);
end

relativeNormalizedW = zeros(size(relativeW));
for ii = 1:length(R2)
    for jj=1:length(R2)
        relativeNormalizedW(ii,jj) = relativeW(ii,jj) * sqrt(varCov(jj,jj))/stdDev(ii);
    end
    idioW(ii,ii) = idioW(ii,ii) / stdDev(ii);
end
if strcmp(decompositionMethod,'cholesky')
    %use cholesky decomposition
    transformMatrix = chol(corrMatrix,'lower');
elseif strcmp(decompositionMethod,'SVD')
    %use singular value decomposition
    [U,S,V] = svd(corrMatrix);
    transformMatrix = (U * S^0.5);
end

relativeWSyntheticFactors = relativeNormalizedW;

%% We now proceed to our simulation
%m = 1000000; %number of simulations
n = length(R2); %Number of obligors in ptf
k = length(indexName); %Number of indices

rng(2262); %seed for synthetic factors
syntheticFactors = (transformMatrix*randn(nbSim,k)')'; %Matrix of simulated N(0,1) for synthetic factors

rng(20119) %seed for noises
noises = randn(nbSim,n); %Matrix of simulated N(0,1) for idiosyncratic noise

rng('default') %seed for recovery rate
recovery = rand(nbSim,n); %Matrix of U(0,1) that simulates probability used for recovery calculation
 %parameters for recovery calculation

%% Compute Returns using simulated values and weights
weightSyntheticAndIdiosyncratic = [relativeWSyntheticFactors idioW];

returns = (weightSyntheticAndIdiosyncratic * [syntheticFactors noises]')';

%Define qualitative ratings associated with simulated returns
horizonRatings = DefineHorizonRating_v2(returns,thresholds,ratings,ptfRatings);

values = transition(2:end,:)/100*repmat(horizonValues',6,1)';
mhv = sum(values(:,1));
%Calculate alpha and beta parameters for the Beta distribution based on
%input mean reacovery rate and standard deviation of recovery rates.
alpha = (meanRR^2*((1-meanRR)/stdRR^2)) - meanRR;
beta = alpha * (1-meanRR)/meanRR;

HorizonPricing = zeros(size(returns,1),size(returns,2)+1);
for ii = 1:nbSim
    for jj = 1:n
        if strcmp(horizonRatings{ii,jj},'Default')
           
            HorizonPricing(ii,jj) = betainv(recovery(ii,jj),alpha,beta)*recoveryClaim;
        else 
            HorizonPricing(ii,jj) = horizonValues(find(strcmp([ratings],horizonRatings{ii,jj})));
        end
    end
    HorizonPricing(ii,end) = sum(HorizonPricing(ii,1:end-1));
end


%% Retreive CreditVaR and plot results
histfit(HorizonPricing(:,end)/1000000);
%mhv = mean(HorizonPricing(:,end));
title(['Distribution of porfolio value at 1 year risk horizon ' num2str(nbSim) ' simulations']);
xlabel('Portfolio value at 1 year risk horizon in million USD');
ylabel('Frequency');

clc;
confidenceLevel = [0.95;0.97;0.99;0.9975;0.9999];
fprintf('*****************************STATS***************************\n');
for i = 1:length(confidenceLevel)
    cVaR = mhv-prctile(HorizonPricing(:,end),(1-confidenceLevel(i))*100);
    %fprintf('************************STATS**********************\n');
    %fprintf('%.2f percent Portfolio value is %.2f $\n',(1-confidenceLevel(i))*100,prctile(HorizonPricing(:,end),(1-confidenceLevel(i))*100));
    fprintf('%.2f%% C-VaR on 1 year horizon: %.2f$ or %.2f percent\n',confidenceLevel(i)*100,cVaR,cVaR/mhv*100);
    %fprintf('****************************************************\n');
end
fprintf('**************************************************************\n');