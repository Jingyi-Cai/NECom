function model_new=make_sub_model(model,sel_rxn_ind,sel_met_ind,community_flag)

%{ 
Function: extract a part from an existing model as a new model

Input: 
model:             Stoichiometric model data struct
sel_rxn_ind:       indices of reactions to be extracted
sel_met_ind        indices of metabolites to be extracted

Output:
model_new          extracted model

Jingyi Cai 2018.2
%}


if nargin < 2
sel_rxn_ind=1:length(model.rxns); % default: all reactions
end


if nargin < 3
sel_met_ind=1:length(model.mets); % default: all metabolites
end


if nargin < 4
community_flag=0; % default: the model is not a community model
end
model_new.rxns=model.rxns(sel_rxn_ind);
model_new.mets=model.mets(sel_met_ind);
model_new.lb=model.lb(sel_rxn_ind);
model_new.ub=model.ub(sel_rxn_ind);
model_new.S=model.S(sel_met_ind,sel_rxn_ind);
model_new.c=model.c(sel_rxn_ind);
if isfield(model,'rev')
    model_new.rev=model.rev(sel_rxn_ind);
end
if community_flag==1
model_new.spBm=model.spBm;
model_new.rxnSps=model.rxnSps(model.rxnSps>=1);
model_new.metSps=model.metSps(model.metSps>=1);
    if isfield(model,'sps')
    model_new.sps=model.sps;
    end

    if isfield(model,'spNames')
    model_new.spNames=model.spNames;
    end
    if isfield(model,'EXsp')
    model_new.EXsp=model.EXsp;
    end

    if isfield(model,'UTsp')
    model_new.UTsp=model.UTsp;
    end    
end


if isfield(model,'rxnNames')& length(model.rxnNames)==length(model.rxns)
     if length(model.rxnNames)==length(model.rxns)    
    model_new.rxnNames=model.rxnNames(sel_rxn_ind);
     end
end

if isfield(model,'metNames')& length(model.metNames)==length(model.mets)
    if length(model.metNames)==length(model.mets) 
    model_new.metNames=model.metNames(sel_met_ind);
    end
end
if isfield(model,'metFormulas')& length(model.metFormulas)==length(model.mets)
    if length(model.metFormulas)==length(model.mets)     
    model_new.metFormulas=model.metFormulas(sel_met_ind);
    end
end

end





