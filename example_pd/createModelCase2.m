
initCobraToolbox
rxnSps=[];metSps=[];sps={'sp1','sp2'};EXsp=[];
themodel1 = createModel; % initialize Cobra toolbox and cobra model 
% adding reactions for sp1
themodel1 = addReaction(themodel1,'EX_S_sp1','reactionFormula','S1[e]_sp1  <=> ','lowerBound',-10,'upperBound',0); % limiting substrate uptake from medium
themodel1 = addReaction(themodel1,'EX_A_sp1','reactionFormula','A[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel1 = addReaction(themodel1,'EX_B_sp1','reactionFormula','B[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel1 = addReaction(themodel1,'EX_P_sp1','reactionFormula','P1[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel1 = addReaction(themodel1,'EX_h_sp1','reactionFormula','h[e]_sp1  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'EX_pi_sp1','reactionFormula','pi[e]_sp1  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'EX_h2o_sp1','reactionFormula','h2o[e]_sp1  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_S_sp1','reactionFormula','S1[e]_sp1  <=> S1[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_A_sp1','reactionFormula','A[e]_sp1  <=> A[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_B_sp1','reactionFormula','B[e]_sp1  <=> B[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_h_sp1','reactionFormula','h[e]_sp1  <=> h[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_P_sp1','reactionFormula','P1[e]_sp1  <=> P1[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_pi_sp1','reactionFormula','pi[e]_sp1  <=> pi[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'T_h2o_sp1','reactionFormula','h2o[e]_sp1  <=> h2o[c]_sp1 ','lowerBound',-Inf,'upperBound',Inf);
themodel1 = addReaction(themodel1,'R_1_sp1','reactionFormula','S1[c]_sp1 + 3 h[c]_sp1 + 3 adp[c]_sp1  -> P1[c]_sp1 + 3 atp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel1 = addReaction(themodel1,'R_2_sp1','reactionFormula','P1[c]_sp1 + atp[c]_sp1  -> B[c]_sp1 + h[c]_sp1 + adp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel1 = addReaction(themodel1,'R_3_sp1','reactionFormula','P1[c]_sp1 + 3 atp[c]_sp1  -> A[c]_sp1 + 3 h[c]_sp1 + 3 adp[c]_sp1 ','lowerBound',0,'upperBound',0);
themodel1 = addReaction(themodel1,'R_4_sp1','reactionFormula','h2o[c]_sp1 + atp[c]_sp1  -> h[c]_sp1 + pi[c]_sp1 + adp[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel1 = addReaction(themodel1,'OBJ_sp1','reactionFormula','3 A[c]_sp1 + 3 B[c]_sp1 + 5 atp[c]_sp1  -> 5 h[c]_sp1 + 5 adp[c]_sp1 + biomass[c]_sp1 ','lowerBound',0,'upperBound',Inf);
themodel1 = addReaction(themodel1,'Biomass_sp1','reactionFormula','biomass[c]_sp1  -> ','lowerBound',0,'upperBound',Inf,'objectiveCoef', 1);
len1=length(themodel1.rxns);len1m=length(themodel1.mets);
rxnSps=[rxnSps;ones(len1,1)];metSps=[metSps;ones(len1m,1)];
%FBA(themodel1)

themodel2 = createModel; % initialize Cobra toolbox and cobra model 
% adding reactions for sp2
themodel2 = addReaction(themodel2,'EX_S_sp2','reactionFormula','S2[e]_sp2  <=> ','lowerBound',-10,'upperBound',0); % limiting substrate uptake from medium
themodel2 = addReaction(themodel2,'EX_A_sp2','reactionFormula','A[e]_sp2  -> ','lowerBound',0,'upperBound',0); 
themodel2 = addReaction(themodel2,'EX_C_sp2','reactionFormula','C[e]_sp2  -> ','lowerBound',0,'upperBound',0);
themodel2 = addReaction(themodel2,'EX_P_sp2','reactionFormula','P[e]_sp2  -> ','lowerBound',0,'upperBound',0);
themodel2 = addReaction(themodel2,'EX_h_sp2','reactionFormula','h[e]_sp2  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'EX_pi_sp2','reactionFormula','pi[e]_sp2  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'EX_h2o_sp2','reactionFormula','h2o[e]_sp2  <=> ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_S_sp2','reactionFormula','S2[e]_sp2  <=> S2[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_A_sp2','reactionFormula','A[e]_sp2  <=> A[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_C_sp2','reactionFormula','C[e]_sp2  <=> C[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_h_sp2','reactionFormula','h[e]_sp2  <=> h[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_P_sp2','reactionFormula','P2[e]_sp2  <=> P2[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_pi_sp2','reactionFormula','pi[e]_sp2  <=> pi[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'T_h2o_sp2','reactionFormula','h2o[e]_sp2  <=> h2o[c]_sp2 ','lowerBound',-Inf,'upperBound',Inf);
themodel2 = addReaction(themodel2,'R_1_sp2','reactionFormula','S2[c]_sp2 + 3 h[c]_sp2 + 3 adp[c]_sp2  -> P2[c]_sp2 + 3 atp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel2 = addReaction(themodel2,'R_2_sp2','reactionFormula','P2[c]_sp2 +   atp[c]_sp2  -> C[c]_sp2 +    h[c]_sp2 +    adp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel2 = addReaction(themodel2,'R_3_sp2','reactionFormula','P2[c]_sp2 + 3 atp[c]_sp2  -> A[c]_sp2 +  3 h[c]_sp2 +  3 adp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel2 = addReaction(themodel2,'R_4_sp2','reactionFormula','h2o[c]_sp2 + atp[c]_sp2  -> h[c]_sp2 + pi[c]_sp2 + adp[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel2 = addReaction(themodel2,'OBJ_sp2','reactionFormula','6 A[c]_sp2 + 6 C[c]_sp2 + 5 atp[c]_sp2  -> 5 h[c]_sp2 + 5 adp[c]_sp2 + biomass[c]_sp2 ','lowerBound',0,'upperBound',Inf);
themodel2 = addReaction(themodel2,'Biomass_sp2','reactionFormula','biomass[c]_sp2  -> ','lowerBound',0,'upperBound',Inf,'objectiveCoef', 1);
len2=length(themodel2.rxns);len2m=length(themodel2.mets);
rxnSps=[rxnSps;2*ones(len2,1)];metSps=[metSps;2*ones(len2m,1)];
%FBA2(themodel2)

%% merge two individual models and setting community exchange rxns-----------for NECom and OptCom (these algorithms has built-in community compartment configuration)
themodel = mergeTwoModels(themodel1,themodel2);
% adding inter-species exchange reactions for sp1
themodel = addReaction(themodel,'EX_A_2_sp1','reactionFormula','A[e]_sp1  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_B_2_sp1','reactionFormula','B[e]_sp1  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_P_2_sp1','reactionFormula','P1[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_S_2_sp1','reactionFormula','S1[e]_sp1  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_h2o_2_sp1','reactionFormula','h2o[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_h_2_sp1','reactionFormula','h[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_pi_2_sp1','reactionFormula','pi[e]_sp1  -> ','lowerBound',0,'upperBound',Inf);
len1ex=length(themodel.rxns);len1mex=length(themodel.mets);EXsp=[EXsp;(len2+1:len1ex)'];
rxnSps1=[rxnSps;ones(len1ex-len2-len1,1)];metSps1=[metSps;ones(len1mex-len2m-len1m,1)];
% adding inter-species exchange reactions for sp2
themodel = addReaction(themodel,'EX_A_2_sp2','reactionFormula','A[e]_sp2  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_C_2_sp1','reactionFormula','C[e]_sp1  <=> ','lowerBound',-5,'upperBound',Inf);
themodel = addReaction(themodel,'EX_P_2_sp2','reactionFormula','P2[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_S_2_sp2','reactionFormula','S2[e]_sp2  -> ','lowerBound',0,'upperBound',0);
themodel = addReaction(themodel,'EX_h2o_2_sp2','reactionFormula','h2o[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_h_2_sp2','reactionFormula','h[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
themodel = addReaction(themodel,'EX_pi_2_sp2','reactionFormula','pi[e]_sp2  -> ','lowerBound',0,'upperBound',Inf);
len2ex=length(themodel.rxns);len2mex=length(themodel.mets);EXsp=[EXsp;(len1ex+1:len2ex)'];
rxnSps1=[rxnSps1;2*ones(len2ex-len1ex,1)];metSps1=[metSps1;2*ones(len2mex-len1mex,1)];

%-------------------------------------------------------------------------
[themodel.rxnSps,themodel.metSps,themodel.EXsp,themodel.sps]=deal(rxnSps1,metSps1,EXsp,sps);
themodel.c=zeros(length(themodel.rxns),1);
themodel.c(findRxnIDs(themodel,{'OBJ_sp1','OBJ_sp2'}))=1;themodel.spBm=find(themodel.c);
%% merge two individual models and setting community exchange rxns-----------for Joint-FBA (need to configure the community compartment manually)

themodel0 = mergeTwoModels(themodel1,themodel2);

% adding inter-species exchange reactions for sp1
themodel0 = addReaction(themodel0,'EX_A_2_sp1','reactionFormula','A[e]_sp1  <=> A[e]_com','lowerBound',-5,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_B_2_sp1','reactionFormula','B[e]_sp1  <=> B[e]_com','lowerBound',-5,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_P_2_sp1','reactionFormula','P1[e]_sp1  -> P1[e]_com','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_S_2_sp1','reactionFormula','S1[e]_sp1  -> S1[e]_com','lowerBound',0,'upperBound',0);
themodel0 = addReaction(themodel0,'EX_h2o_2_sp1','reactionFormula','h2o[e]_sp1  -> h2o[e]_com','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_h_2_sp1','reactionFormula','h[e]_sp1  -> h[e]_com','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_pi_2_sp1','reactionFormula','pi[e]_sp1  -> pi[e]_com','lowerBound',0,'upperBound',Inf);
len1ex2=length(themodel0.rxns);len1mex2=length(themodel0.mets);EXsp=[EXsp;(len2+1:len1ex2)'];
rxnSps2=[rxnSps;ones(len1ex2-len2-len1,1)];metSps2=[metSps;zeros(len1mex2-len2m-len1m,1)];
% adding inter-species exchange reactions for sp2
themodel0 = addReaction(themodel0,'EX_A_2_sp2','reactionFormula','A[e]_sp2  <=> A[e]_com','lowerBound',-5,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_C_2_sp2','reactionFormula','C[e]_sp2  <=> C[e]_com','lowerBound',-5,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_P_2_sp2','reactionFormula','P2[e]_sp2  -> P2[e]_com','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_S_2_sp2','reactionFormula','S2[e]_sp2  -> S2[e]_com','lowerBound',0,'upperBound',0);
themodel0 = addReaction(themodel0,'EX_h2o_2_sp2','reactionFormula','h2o[e]_sp2  -> h2o[e]_com','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_h_2_sp2','reactionFormula','h[e]_sp2  -> h[e]_com','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_pi_2_sp2','reactionFormula','pi[e]_sp2  -> pi[e]_com','lowerBound',0,'upperBound',Inf);
len2ex2=length(themodel0.rxns);len2mex2=length(themodel0.mets);EXsp=[EXsp;(len1ex2+1:len2ex2)'];
rxnSps2=[rxnSps2;2*ones(len2ex2-len1ex2,1)];metSps2=[metSps2;zeros(len2mex2-len1mex2,1)];

% adding community level exchange
themodel0 = addReaction(themodel0,'EX_A_com','reactionFormula','A[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_B_com','reactionFormula','B[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_C_com','reactionFormula','C[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_P1_com','reactionFormula','P1[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_P2_com','reactionFormula','P2[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_S1_com','reactionFormula','S1[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_S2_com','reactionFormula','S2[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_h2o_com','reactionFormula','h2o[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_h_com','reactionFormula','h[e]_com  <=> ','lowerBound',0,'upperBound',Inf);
themodel0 = addReaction(themodel0,'EX_pi_com','reactionFormula','pi[e]_com  <=> ','lowerBound',0,'upperBound',Inf);

len0ex=length(themodel0.rxns);len0mex=length(themodel0.mets);
rxnSps2=[rxnSps2;zeros(len0ex-len2ex2,1)];metSps2=[metSps2;zeros(len0mex-len2mex2,1)];
[themodel0.rxnSps,themodel0.metSps,themodel0.EXsp,themodel0.sps]=deal(rxnSps2,metSps2,EXsp,sps);
themodel0.c=zeros(length(themodel0.rxns),1);
themodel0.c(findRxnIDs(themodel0,{'OBJ_sp1','OBJ_sp2'}))=1;themodel0.spBm=find(themodel0.c);

















