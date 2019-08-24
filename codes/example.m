
% demonstrative example for NECom and FShaP simulations
% Jingyi Cai 2019-8
clear;clc;
% load the community model, which consists of two individual models,
load themodel.mat themodel
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

% the metabolites may only be exchanged between sp1 and sp2 are A and B(other metabolite 
% such as proton can be exchanged between organisms and environment, so not considered here)
% let's label A and B a metabolite number 1,2 repectively
% the associated crossfeeding reactions are EX_A_2_sp1 EX_B_2_sp1 EX_A_2_sp2 EX_B_2_sp2
% they are given an order 
parameters.cf_order=[1;2;3;4];
% find their reaction indices in the community model:
parameters.sub_indExSpi=[41;42;48;49];
% amount of crossfeeding reactions:
parameters.numSub_ExRxn=length(parameters.sub_indExSpi);
% since they are crossfeeding reactions for A,B,A,B, so their metabolite number are:  
parameters.spsInds=[1 2 1 2];
% and they belong to sp1,sp1,sp2, sp2, respectively:
parameters.sub_exSpi=[1;1;2;2];
% 
% crossfeeding reactions that can affect the availability of substrate
% uptaken by other species: EX_A_2_sp2 EX_B_2_sp2 EX_A_2_sp1 EX_B_2_sp1 and
% their indices:
parameters.other_Ex_all=[48;49;41;42];
% they belong to sp2,sp2,sp1,sp1 respectively:
parameters.other_sp_ind=[2;2;1;1];
% mapping to the reaction they may affect(refering to cf_order):
parameters.order_other=[1;2;3;4];

%---------------Final state predictions----------------------------
% run Joint-FBA simulation
solfba1=fba(themodel); % Flux Balance Analysis
x1=solfba1.x;
%  joint-FBA only supports abundance ratio of 1:1, in order to compare 
%  different methods in paralell we set Relative abundance for sp1 and sp2:
X=[0.5,0.5]; % they are 1:1
% run OptCom simulation
[x2,exitflag2,index2,Lmodel2,fluxes2,info2]=OptCom(themodel,parameters,X);
% run NECom simulation
[x3,exitflag3,index3,Lmodel3,fluxes3,info3]=NECOM(themodel,parameters,X);
Indrxns=find(themodel.rxnSps>=1);
allfluxes=[x1(Indrxns),x2(Indrxns),x3(Indrxns)];
% if MS Office Excel install we can write it into the excel file. 
% need to uncheck all COM Add-in in Excel options to make it works
try xlswrite('pd.xls',allfluxes,'reactions','H2');
catch
end
%---------------FShaP ----------------------------------------------

% prepare the state TransitionMap, the range of crossfeeding of A and B 
% will be set the maximal
% crossfeeding values are custom-setup, but in this example we set them as the
% solution values of Joint-FBA(optCom has the same values), in order to
% test if Joint-FBA(and OptCom) will predict Nash�?quilibria

maxrates(1)=-x2(parameters.sub_indExSpi(1)); % the absolute flux of EX_A_2_sp1
maxrates(2)=-x2(parameters.sub_indExSpi(4)); % the absolute flux of EX_B_2_sp2
num_ins=[5,5];% there are total 5X5 states on the map
EX_Rflux{1,1}=linspace(0,maxrates(1),num_ins(1)); %  crossfeeding rates of A, or X axis ticks of the map
EX_Rflux{1,2}=linspace(0,maxrates(2),num_ins(2)); %  crossfeeding rates of B, or Y axis ticks of the map
% controling reactions are EX_A_2_sp1 EX_B_2_sp2 for uptake of A, B respectively  
% by controling these two uptake reactions we can capture state for StateTransitionMap
ctrlinds=parameters.sub_indExSpi([3,2]); 
% generate payoff matrices(FShaP) for the state transition map
theresults=StateTransitionMap(themodel,parameters,ctrlinds,EX_Rflux);

% if MS Office Excel install we can write it into the excel file. 
% need to uncheck all COM Add-in in Excel options to make it works
theVarRxns=themodel.rxns(ctrlinds);
EX_Rflux=[EX_Rflux{1,1};EX_Rflux{1,2}];
try printMatrixes(theresults,theVarRxns,EX_Rflux,[4,4],'pd.xls')
catch
end













