function [horizonRating] = DefineHorizonRating_v2( returns,thresholds,ratings,ptfRatings )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns a vector of rating for each return scenario
%   Input:
%       returns: matrix of horizon returns for an obligor NxM N:scenarios,
%                M:obligors
%       thresholds: matrix of Zscores for given rating system
%       ratings: cell array containing letters ratings i.e: AAA,AA, etc.
%       ptfRatings: ratings of obligors in portfolio, has to be in order of
%                   their returns
%   Output:
%       rating: Vector of mapped rating corresponding to each return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preallocation
[m,n] = size(returns);
[r,c] = size(thresholds);
horizonRating = cellstr(char(zeros(size(returns)))); %Allocate space for string ratings

%%
for ii = 1:length(ptfRatings) %Parse all obligors in ptf
    % Find the right threshold row for this obligor
    obligorThresholds = thresholds(find(strcmp([ratings],ptfRatings{ii})),:);
    obligorReturns = returns(:,ii);
    
    %Treat default seperatly
    horizonRating(find(obligorReturns<=obligorThresholds(end)),ii) = ratings(end);
    
    %Treat all other thresholds
    for jj = 2:length(obligorThresholds)
        firstFilter = find(obligorReturns<obligorThresholds(jj-1));
        secondFilter = find(obligorReturns>=obligorThresholds(jj));
        
        horizonRating(find(obligorReturns<obligorThresholds(jj-1) & obligorReturns>=obligorThresholds(jj)),ii) = ...
            ratings(jj-1);
    end   
end
% %on boucle sur tous les obligors
% for ii = 1:n
%   %On doit traiter le premier et dernier threshold séparément
%   numRating(find(returns(:,ii)> thresholds(1,ii)),ii) = 1; %Premier threshold
%   numRating(find(returns(:,ii) < thresholds(end,ii)),ii) = length(ratings);
%    
%     for jj = 1:r-1
%       numRating(find((returns(:,ii)< thresholds(jj,ii)) & (returns(:,ii) >...
%       thresholds(jj+1,ii))),ii) = jj+1;
%     end
% end
% 
% %% Mapping ratings from numerical to actual string ratings
% for ll = n:-1:1 %parse all obligors
%     for kk = 1:length(ratings) %parse all ratings
%         rating(find(numRating(:,ll)==kk),ll) = ratings(kk);
%     end
% end

end

