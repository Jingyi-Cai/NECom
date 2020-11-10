% 
% demonstrative example for NECom and simulations
% This example is requested by a reviewer. 
% Jingyi Cai 2020-10
clear;clc;
thesolver='gurobi'; % choose from 'baron','bonmin' and 'gurobi', bonmin solver requires installation of OPTI toolbox on Windows 
% load the community model, which consists of two individual models,
%load ../themodel.mat themodel

% or directly create the community model using CobraToolbox

createModelCase2



% spreadsheet version of this model presents in 'Reactions',and
% 'Metabolites' page in pd.xls
%{
Data struct of the model:
(Assume the community model has I metabolites and J reactions, belonging to
N species. )

       rxns:       J by 1 string cell reactions
       mets:       I by 1 string cell metabolites
       S:          I by J Stoichiometry matrix of the model
       c:          J by 1 objective 0-1 vector 
       lb:         J by 1 real number vector for lower bounds of reactions 
       ub:         J by 1 real number vector for upper bounds of reactions 
       rxnSps:     J by 1 vector indicating to which species a reaction belong
                    (modelCom.rxns{j} belongs to species #modelCom.rxnSps(j))
       metSps:     I by 1 vector indicating to which species a metabolite belong
                    (modelCom.mets{j} belongs to species #modelCom.rxnSps(j),
                    0 for community metabolites)
       spBm:       1 by N vector indicating the ID of biomass reaction for each species
                    (species k has biomass reaction modelCom.rxns{modelCom.spBm(k)})
       EXsp:       Ids for exchange reactions between species and
                   community [a_ij] is the id for the export of i-th
%}

% now configure parameters for simulation

% the metabolites may only be exchanged between sp1 and sp2 is A(other metabolite 
% such as proton can be exchanged between organisms and environment, so not considered here)
% let's label A a metabolite number 1
% the associated crossfeeding reactions are EX_A_2_sp1  EX_A_2_sp2 
% they are given an order 
parameters.cf_order=[1;2];
% find their reaction indices in the community model:
parameters.sub_indExSpi=[41;48];
% number of crossfeeding reactions:
parameters.numSub_ExRxn=length(parameters.sub_indExSpi);
% since they are crossfeeding reactions for A,A so their metabolite number are:  
parameters.spsInds=[1 1];
% and they belong to sp1,sp2 respectively:
parameters.sub_exSpi=[1;2];
% 
% crossfeeding reactions that can affect the availability of substrate
% uptaken by other species: EX_A_2_sp2  EX_A_2_sp1  and
% their indices:
parameters.other_Ex_all=[48;41];
% they belong to sp2,sp2,sp1,sp1 respectively:
parameters.other_sp_ind=[2;1];
% mapping to the reaction they may affect(refering to cf_order):
parameters.order_other=[1;2];

%---------------Final state predictions----------------------------
% run Joint-FBA simulation
[solfba1,Lmodel1,index1]=fba(themodel0); % Flux Balance Analysis
x1=solfba1.x;
%  joint-FBA only supports abundance ratio of 1:1, in order to compare 
%  different methods in paralell we set Relative abundance for sp1 and sp2:
X=[0.5,0.5]; % they are 1:1
% run OptCom simulation
[x2,exitflag2,index2,Lmodel2,fluxes2,info2]=OptCom(themodel,parameters,X,thesolver);
% run NECom simulation
[x3,exitflag3,index3,Lmodel3,fluxes3,info3]=NECOM(themodel,parameters,X,thesolver);
Indrxns=find(themodel.rxnSps>=1);
allfluxes=[x1(Indrxns),x2(Indrxns),x3(Indrxns)];
% if MS Office Excel install we can write it into the excel file. 
% need to uncheck all COM Add-in in Excel options to make it works
writeCbModel(themodel,'xls','modeltext.xls')
try xlwrite('modeltext.xls',allfluxes,'Reaction List','H2');xlwrite('modeltext.xls',{'JointFBA','OptCom','NECom'},'Reaction List','H1');
    fprintf('fluxes predicted by jointFBA, OptCom and NEcom has been saved in modeltext.xls or modeltext.csv\n')
catch
    fprintf('fluxes predicted by jointFBA, OptCom and NEcom can be found in allfluxes matrix\n')
end


%{
fprintf('Test if the flux predicted by JointFBA is feasible to NECom...');
externalconstr.lbindices=1:length(themodel.rxns);
externalconstr.ubindices=1:length(themodel.rxns);
externalconstr.lb=x1(Indrxns);
externalconstr.ub=x1(Indrxns);
[xtest1,exitflagtest1,indextest1,Lmodeltest1,fluxestest1,infotest1]=NECOM(themodel,parameters,X,thesolver,externalconstr,'FBATest');
fprintf(['The flux predicted by JointFBA is ',infotest1, ' to NECom\n']);

fprintf('Test if the flux predicted by OptCom is feasible to NECom...');
externalconstr.lbindices=1:length(themodel.rxns);
externalconstr.ubindices=1:length(themodel.rxns);
externalconstr.lb=x2(Indrxns);
externalconstr.ub=x2(Indrxns);
[xtest2,exitflagtest2,indextest2,Lmodeltest2,fluxestest2,infotest2]=NECOM(themodel,parameters,X,thesolver,externalconstr,'OptComTest');
fprintf(['The flux predicted by OptCom is ',infotest2, ' to NECom\n']);

%------------To determine the minimal subset of constraints (including bounds) responsible for
%the infeasibility the model should be exported to python environment,
%where we can call Model.computellS() to determine the minimal subset of
%infeasible constraints and bounds
fprintf('please open jupyter notebook and run checkInfeasibility.ipynb to find the minimal subset of constraints and bounds that JointFBA solution fails to satisfy\n')
%------------Further Test-----------------------------------------------
%}







































