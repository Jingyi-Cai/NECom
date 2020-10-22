
initCobraToolbox
rxnSps=[];metSps=[];sps={'sp1','sp2'};EXsp=[];
themodel = createModel; % initialize Cobra toolbox and cobra model 
% adding reactions for sp1
themodel = addReaction(themodel,'EX_S_sp1','reactionFormula','S[e]_sp1  <=> ','lowerBound',-10,'upperBound',0); % limiting substrate uptake from medium
themodel = addReaction(themodel,'EX_A_sp1','reactionFormula','A[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_B_sp1','reactionFormula','B[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_P_sp1','reactionFormula','P[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_h_sp1','reactionFormula','h[e]_sp1  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'EX_pi_sp1','reactionFormula','pi[e]_sp1  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'EX_h2o_sp1','reactionFormula','h2o[e]_sp1  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_S_sp1','reactionFormula','S[e]_sp1  <=> S[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_A_sp1','reactionFormula','A[e]_sp1  <=> A[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_B_sp1','reactionFormula','B[e]_sp1  <=> B[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_h_sp1','reactionFormula','h[e]_sp1  <=> h[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_P_sp1','reactionFormula','P[e]_sp1  <=> P[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_pi_sp1','reactionFormula','pi[e]_sp1  <=> pi[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_h2o_sp1','reactionFormula','h2o[e]_sp1  <=> h2o[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'R_1_sp1','reactionFormula','S[c]_sp1 + 2 h[c]_sp1 + 2 adp[c]_sp1  -> P[c]_sp1 + 2 atp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'R_2_sp1','reactionFormula','P[c]_sp1 + atp[c]_sp1  -> B[c]_sp1 + h[c]_sp1 + adp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'R_3_sp1','reactionFormula','P[c]_sp1 + 3 atp[c]_sp1  -> A[c]_sp1 + 3 h[c]_sp1 + 3 adp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'R_4_sp1','reactionFormula','h2o[c]_sp1 + atp[c]_sp1  -> h[c]_sp1 + pi[c]_sp1 + adp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'OBJ_sp1','reactionFormula','3 A[c]_sp1 + 3 B[c]_sp1 + 5 atp[c]_sp1  -> 5 h[c]_sp1 + 5 adp[c]_sp1 + biomass[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'Biomass_sp1','reactionFormula','biomass[c]_sp1  -> ','lowerBound',0,'upperBound',Inf,'objectiveCoef', 1);
len1=length(themodel.rxns);spBm(1,1)=len1;len1m=length(themodel.mets);
rxnSps=[rxnSps;ones(len1,1)];metSps=[metSps;ones(len1m,1)];


% adding reactions for sp2
themodel = addReaction(themodel,'EX_S_sp2','reactionFormula','S[e]_sp2  <=> ','lowerBound',-10,'upperBound',0); % limiting substrate uptake from medium
themodel = addReaction(themodel,'EX_A_sp2','reactionFormula','A[e]_sp2  -> ','lowerBound',0,'upperBound',0); 
themodel = addReaction(themodel,'EX_B_sp2','reactionFormula','B[e]_sp2  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_P_sp2','reactionFormula','P[e]_sp2  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_h_sp2','reactionFormula','h[e]_sp2  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'EX_pi_sp2','reactionFormula','pi[e]_sp2  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'EX_h2o_sp2','reactionFormula','h2o[e]_sp2  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_S_sp2','reactionFormula','S[e]_sp2  <=> S[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_A_sp2','reactionFormula','A[e]_sp2  <=> A[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_B_sp2','reactionFormula','B[e]_sp2  <=> B[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_h_sp2','reactionFormula','h[e]_sp2  <=> h[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_P_sp2','reactionFormula','P[e]_sp2  <=> P[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_pi_sp2','reactionFormula','pi[e]_sp2  <=> pi[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'T_h2o_sp2','reactionFormula','h2o[e]_sp2  <=> h2o[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel = addReaction(themodel,'R_1_sp2','reactionFormula','S[c]_sp2 + 2 h[c]_sp2 + 2 adp[c]_sp2  -> P[c]_sp2 + 2 atp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'R_2_sp2','reactionFormula','P[c]_sp2 + 3 atp[c]_sp2  -> B[c]_sp2 + 3 h[c]_sp2 + 3 adp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'R_3_sp2','reactionFormula','P[c]_sp2 + atp[c]_sp2  -> A[c]_sp2 + h[c]_sp2 + adp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'R_4_sp2','reactionFormula','h2o[c]_sp2 + atp[c]_sp2  -> h[c]_sp2 + pi[c]_sp2 + adp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'OBJ_sp2','reactionFormula','3 A[c]_sp2 + 3 B[c]_sp2 + 5 atp[c]_sp2  -> 5 h[c]_sp2 + 5 adp[c]_sp2 + biomass[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'Biomass_sp2','reactionFormula','biomass[c]_sp2  -> ','lowerBound',0,'upperBound',Inf,'objectiveCoef', 1);
len2=length(themodel.rxns);spBm(2,1)=len2;len2m=length(themodel.mets);
rxnSps=[rxnSps;2*ones(len2-len1,1)];metSps=[metSps;2*ones(len2m-len1m,1)];


% adding inter-species exchange reactions for sp1
themodel = addReaction(themodel,'EX_A_2_sp1','reactionFormula','A[e]_sp1  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_B_2_sp1','reactionFormula','B[e]_sp1  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_P_2_sp1','reactionFormula','P[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_S_2_sp1','reactionFormula','S[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_h2o_2_sp1','reactionFormula','h2o[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_h_2_sp1','reactionFormula','h[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_pi_2_sp1','reactionFormula','pi[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
len1ex=length(themodel.rxns);len1mex=length(themodel.mets);EXsp=[EXsp,(len2+1:len1ex)'];
rxnSps=[rxnSps;ones(len1ex-len2,1)];metSps=[metSps;ones(len1mex-len2m,1)];
% adding inter-species exchange reactions for sp2
themodel = addReaction(themodel,'EX_A_2_sp2','reactionFormula','A[e]_sp2  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_B_2_sp2','reactionFormula','B[e]_sp2  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_P_2_sp2','reactionFormula','P[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_S_2_sp2','reactionFormula','S[e]_sp2  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_h2o_2_sp2','reactionFormula','h2o[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_h_2_sp2','reactionFormula','h[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_pi_2_sp2','reactionFormula','pi[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
len2ex=length(themodel.rxns);len2mex=length(themodel.mets);EXsp=[EXsp,(len1ex+1:len2ex)'];
rxnSps=[rxnSps;2*ones(len2ex-len1ex,1)];metSps=[metSps;2*ones(len2mex-len1mex,1)];

%-------------------------------------------------------------------------
[themodel.rxnSps,themodel.metSps,themodel.EXsp,themodel.sps,themodel.spBm]=deal(rxnSps,metSps,EXsp,sps,spBm);















