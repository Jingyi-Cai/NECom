function [thematrix,StrategyFeasibilitys]=FShaP(themodel,parameters,EX_R_flux)

%{
Function: construct payoff matrix 
Input: 
    themodel: microbial community model struct
    parameters: parameters 
    EX_R_flux: crossfeeding fluxes to define current interaction state. 
Output: 
    thematrix: payoff matrix(FSHaP)
    StrategyFeasibilitys: strategy feasibility indicator matrix for the
    payoff matrix, 1 indicate the corresponding strategy is feasible, NaN
    otherwise
Jingyi Cai 2019-03-02
%}

subModels=[];
solutions={};
number_sps=length(themodel.spBm);

%fix the lb and ub of the crossfeeding reactions to the crossfeeding fluxes
themodel_1=themodel;
themodel_1.lb(parameters.sub_indExSpi)=EX_R_flux;
themodel_1.ub(parameters.sub_indExSpi)=EX_R_flux;

% apply current Ex bounds to the modelEx

selExInds=full(sparse(parameters.spsInds,parameters.sub_exSpi,...
parameters.sub_indExSpi,length(unique(parameters.sub_exSpi)),length(unique(parameters.spsInds))));
impact=[];
for k=1:number_sps
    rxnsInds{k,1}=find(themodel.rxnSps==k);
    [~,selRxnInds{k,1}]=ismember(selExInds(:,k),rxnsInds{k,1});
    subModels{k,1}=make_sub_model(themodel_1,find(themodel_1.rxnSps==k),find(themodel_1.metSps==k)); 
    
    [solutions{k,1}]=fba(subModels{k,1},'all');
    lbraw{k,1}=solutions{k,1}.mu_lb(selRxnInds{k,1});
    ubraw{k,1}=solutions{k,1}.mu_ub(selRxnInds{k,1});
end

%generating strategy combo for both species
[MCom1] = flipud(permn([0,1], sum(parameters.sub_exSpi==1)));
[MCom2] = flipud(permn([0,1], sum(parameters.sub_exSpi==2)));

theeps=1e-3;
StrategyFeasibilitys=[];
MComs{1,1}=MCom1;MComs{2,1}=MCom2;
% determine the payoff matrix size
m=size(MComs{1,1},1);
n=length([0,1]);

for k=1:m
    StrategyFeasibility=[];
    for i=1:number_sps
         StrategyFeasibility{i,1}=ones(1,2);

        for p=1:n
            if MComs{i,1}(k,p)
                themodel0{i,1}=subModels{i};
                % then we need to do some test to see if this strategy is
                % available for both species
                % when calling a strategy(secretion /uptake one unit more) but it is not feasible, then NaN will be assigned    

                % check the secretion ability for each strategy
                themodel0{i,1}.lb(selRxnInds{i,1}(p))=solutions{i,1}.x((selRxnInds{i,1}(p)))+theeps;
                themodel0{i,1}.ub(selRxnInds{i,1}(p))=solutions{i,1}.x((selRxnInds{i,1}(p)))+theeps;
                thesol=fba(themodel0{i,1});
                if isempty(thesol.x)
                    StrategyFeasibility{i,1}(1)=NaN;
                end                
                % check the uptake ability
                themodel0{i,1}.lb(selRxnInds{i,1}(p))=solutions{i,1}.x((selRxnInds{i,1}(p)))-theeps;
                themodel0{i,1}.ub(selRxnInds{i,1}(p))=solutions{i,1}.x((selRxnInds{i,1}(p)))-theeps;
                thesol=fba(themodel0{i,1});
                if isempty(thesol.x)
                    StrategyFeasibility{i,1}(2)=NaN;
                end

            end

        end
    end
    StrategyFeasibilitys=[StrategyFeasibilitys,StrategyFeasibility];
end
StrategyFeasibilitys=StrategyFeasibilitys';

%    name1={'11','10','01','00'};
%    name2={'11','10','01','00'};

for k=1:number_sps
    theothersps=setdiff([1,2],k);
    impact{k,1}= [-solutions{k,1}.mu_lb(selRxnInds{k,1}),solutions{k,1}.mu_lb(selRxnInds{theothersps,1})];
    totalimpact{k,1}=MComs{k,1}*impact{k,1};
end




    thematrix=cell(m,n);
for i=1:m
    for j=1:n
        % secretion payoff for species 1
       if any(MComs{1}(i,:)) % if species 1 currently has secretion strategy....
           comp1= totalimpact{1}(i,1)*StrategyFeasibilitys{i,1}(1);  % calculate the payoffs
       else
           comp1=0;
       end
       % uptake payoff for species 1
       if any(MComs{2}(j,:)) % if species 2 has secrete strategy...
           comp2= totalimpact{1}(j,2)*StrategyFeasibilitys{j,1}(2);
       else
           comp2=0;
       end
       % secretion payoff for species 2
       if any(MComs{2}(j,:)) % if species 2 has secrete strategy...
           comp3= totalimpact{2}(j,1)*StrategyFeasibilitys{j,2}(1);
       else
           comp3=0;
       end
       % uptake payoff for species 2
       if any(MComs{1}(i,:)) % if species 1 has secrete strategy...
           comp4= totalimpact{2}(i,2)*StrategyFeasibilitys{i,2}(2); 
       else
           comp4=0;
       end
       thematrix{i,j}= [round(comp1+comp2,2),...
           round(comp3+comp4,2)];
    end
end

    
    
    



