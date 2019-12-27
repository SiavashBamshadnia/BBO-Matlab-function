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

clc;
clear;
close all;

%% Problem Definition

CostFunction=@(x) Sphere(x);        % Cost Function

nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=-10;         % Decision Variables Lower Bound
VarMax= 10;         % Decision Variables Upper Bound

%% BBO Parameters

MaxIt=1000;          % Maximum Number of Iterations

nPop=50;            % Number of Habitats (Population Size)

KeepRate=0.2;                   % Keep Rate
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats

nNew=nPop-nKeep;                % Number of New Habitats

% Migration Rates
mu=linspace(1,0,nPop);          % Emmigration Rates
lambda=1-mu;                    % Immigration Rates

alpha=0.9;

pMutation=0.1;

sigma=0.02*(VarMax-VarMin);

[BestSol,BestCost]=BBO(CostFunction,nVar,VarSize,VarMin,VarMax,MaxIt,nPop,KeepRate,nKeep,nNew,mu,lambda,alpha,pMutation,sigma)

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
