%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA113
% Project Title: Biogeography-Based Optimization (BBO) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% ________________________________________________________________
% Modified by Siavash Bamshadnia (amtr.ir)
% Email: sbamtr@gmail.com
% For more information, please visit the following link:
% https://uk.mathworks.com/matlabcentral/fileexchange/52901-biogeography-based-optimization-bbo
% ________________________________________________________________

%
%% Problem Definition
% CostFunction=Your Cost Function
% nVar=Number of Decision Variables
% VarSize=Decision Variables Matrix Size
% VarMin=Decision Variables Lower Bound
% VarMax=Decision Variables Upper Bound
%% BBO Parameters
% MaxIt=Maximum Number of Iterations
% nPop=Number of Habitats (Population Size)
% KeepRate=Keep Rate
% nKeep=Number of Kept Habitats
% nNew=Number of New Habitats
% mu=Emmigration Rates
% lambda=Immigration Rates
% alpha=Alpha;
% pMutation=Permutation;
% sigma=Sigma
%

function[BestSol,BestCost]=BBO(CostFunction,nVar,VarSize,VarMin,VarMax,MaxIt,nPop,KeepRate,nKeep,nNew,mu,lambda,alpha,pMutation,sigma)
% Empty Habitat
habitat.Position=[];
habitat.Cost=[];

% Create Habitats Array
pop=repmat(habitat,nPop,1);

% Initialize Habitats
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
end

% Sort Population
[~, SortOrder]=sort([pop.Cost]);
pop=pop(SortOrder);

% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop

for it=1:MaxIt
    
    newpop=pop;
    for i=1:nPop
        for k=1:nVar
            % Migration
            if rand<=lambda(i)
                % Emmigration Probabilities
                EP=mu;
                EP(i)=0;
                EP=EP/sum(EP);
                
                % Select Source Habitat
                j=RouletteWheelSelection(EP);
                
                % Migration
                newpop(i).Position(k)=pop(i).Position(k) ...
                    +alpha*(pop(j).Position(k)-pop(i).Position(k));
                
            end
            
            % Mutation
            if rand<=pMutation
                newpop(i).Position(k)=newpop(i).Position(k)+sigma*randn;
            end
        end
        
        % Apply Lower and Upper Bound Limits
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
        % Evaluation
        newpop(i).Cost=CostFunction(newpop(i).Position);
    end
    
    % Sort New Population
    [~, SortOrder]=sort([newpop.Cost]);
    newpop=newpop(SortOrder);
    
    % Select Next Iteration Population
    pop=[pop(1:nKeep)
         newpop(1:nNew)];
     
    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Update Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
end