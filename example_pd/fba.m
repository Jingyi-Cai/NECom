function [thesolution,Lmodel,index]=fba(model,pridua,thesense,solverset)
% Function: Flux Balance Analysis

% Input: 

% model: the model struct at least consists of fields of metabolite list 'mets', reaction list 'rxns', 
% stoichiometric matrix 'S', flux bounds 'lb','ub',and objective vector 'c'
% thesense: either 'max'(default) or 'min'
% solverset: at least contains field 'thesolver', which specify which
% solver is used for solving fba problem, glpk as default. 

% Output:
% thesolution: the solution of FBA
% Lmodel: model struct of FBA
% index: indices of variables in the FBA problem

% jingyi cai 2017-8
% jingyi cai 2018-12 merge primal and dual solutions
if nargin<2
pridua='primal';
end

if nargin<3
    thesense='maximize';
end
if nargin<4
   if exist('glpk')
    solverset.thesolver='glpk';
   elseif exist('Cplex')
       solverset.thesolver='cplex';
   elseif exist('opti')
       solverset.thesolver='opti';
   else
       error('this fba currently only supports glpk,cplex,opti, you may upgrade the code if you have other solver')
   end
end

% make the model bounded;
bigM=1000;
%reactions with unlimited lb 
lbInf1 = model.lb <= -(bigM-1);
model.lb(lbInf1)=-bigM;
%reactions with unlimited ub
ubInf1 = model.ub >=  (bigM-1);
model.ub(ubInf1)= bigM;

% dual objective sense is always opposite to primal problem's
switch thesense
    case 'maximize'
    thesense2='minimize';
    case 'minimize'
    thesense2='maximize';    
end

switch pridua
    case 'primal'
        [thesolution,Lmodel,index]=FBA_primal(model,thesense,solverset);
        thesolution.problemtype='primal';
    case 'dual'
        [thesolution,Lmodel,index]=FBA_dual(model,thesense2,solverset);
        thesolution.problemtype='dual';
    case 'all'
        [thesolutionp,Lmodelp,indexp]=FBA_primal(model,thesense,solverset);
        [thesolutiond,Lmodeld,indexd]=FBA_dual(model,thesense2,solverset);
        thesolution=CombineStructs(thesolutionp,thesolutiond);
        Lmodel.p=Lmodelp;
        Lmodel.d=Lmodeld;
        index.p=indexp;
        index.d=indexd;
        thesolution.problemtype='primal-dual';    
end

end

function [thesolution,Lmodel,index]=FBA_primal(model,thesense,solverset)

% Function: solving Flux Balance Analysis primal problem
% jingyi cai 2017-7

Lmodel = start_Lmodel(); % initial modelling
% number of fluxes:
n=length(model.rxns);
% number of mets:
m=length(model.mets);
% add variables: flux-->v
varnames=char(strcat('v_',model.rxns));
Lmodel=add_variables(Lmodel,model.c(:), model.lb, model.ub,char(ones(1,length(model.lb))*'C'),varnames);
history_length=0;
history_length_v_1=history_length;
current_length_v_1=history_length_v_1+n;
index.var.v1 = 1:n;
% Sv = 0
Lmodel=add_constraints(Lmodel,zeros(m,1), model.S, zeros(m,1)); 
index.con.mb1 = 1:m;
nCon = index.con.mb1(end);
Lmodel.Model.sense=thesense;

sol_lp=solveLinearProblem(Lmodel,solverset);

if isempty(sol_lp.x)
    thesolution.x=[];
    thesolution.f_p=[];
    thesolution.status=sol_lp.status;
else
    thesolution.x=sol_lp.x(1:n);
    thesolution.f_p=sol_lp.objval;
    thesolution.status=sol_lp.status;
end
%end
thesolution.solver=sol_lp.solver;
thesolution.rxns=model.rxns;
thesolution.mets=model.mets;
end

function [solution_d,Lmodel_d,index]=FBA_dual(model,thesense,solverset)
% Function: solving Flux Balance Analysis dual problem


% jingyi cai 2017-7

if nargin<2
    thesense='minimize';
end
if nargin<3
   if exist('glpk')
    solverset.thesolver='glpk';
   elseif exist('Cplex')
       solverset.thesolver='cplex';
   elseif exist('opti')
       solverset.thesolver='opti';
   else
       error('this fba currently only supports glpk,cplex,opti, you may upgrade the code if you have other solver')
   end
end

Lmodel_d = start_Lmodel(); % initial modelling
M=1000;
% number of flux:
n=length(model.rxns);
m=length(model.mets);
%add lambda
Lmodel_d=add_variables(Lmodel_d,zeros(m,1),-M*ones(m,1), M*ones(m,1));
history_length=0;
history_length_lambda_1=history_length;
current_length_lambda_1=history_length_lambda_1+m;
index.var.lam1 = 1 : m;
nVar = index.var.lam1(end);

%add mu_LB
Lmodel_d=add_variables(Lmodel_d,-model.lb, zeros(n,1), M*ones(n,1));
history_length=current_length_lambda_1;
history_length_mu_LB_1=history_length;
current_length_mu_LB_1=history_length_mu_LB_1+n;
index.var.muLP= nVar + 1:nVar + n;
nVar = nVar + n;

%add mu_UB
Lmodel_d=add_variables(Lmodel_d,model.ub,zeros(n,1), M*ones(n,1));
history_length=current_length_mu_LB_1;
history_length_mu_UB_1=history_length;
current_length_mu_UB_1=history_length_mu_UB_1+n;
index.var.muUP= nVar + 1:nVar + n;
nVar = nVar + n;

nCon=0;
%S'lambda_1 + mu_UB_1 - mu_LB_1 = c
Lmodel_d=add_constraints(Lmodel_d,model.c(:), [model.S', -speye(n), speye(n)], model.c(:));
index.con.dv1 = 1 : nCon + n;
nCon = nCon + n;


Lmodel_d.Model.sense=thesense;
sol_lp_d=solveLinearProblem(Lmodel_d,solverset);
if isempty(sol_lp_d.x)
    solution_d.mu_lb=[];
    solution_d.mu_ub=[];
    solution_d.f_d=[];
    solution_d.status=sol_lp_d.status;
else
    solution_d.lambda=sol_lp_d.x(1:m);
    solution_d.mu_lb=sol_lp_d.x(m+1:m+n);
    solution_d.mu_ub=sol_lp_d.x(m+n+1:m+2*n);
    solution_d.f_d=sol_lp_d.objval; 
    solution_d.status=sol_lp_d.status;
end
solution_d.solver=sol_lp_d.solver;
solution_d.rxns=model.rxns;
solution_d.mets=model.mets;


end





function sol_lp=solveLinearProblem(Lmodel,solverset)

% Function: Solve Linear Planing Problem

% Input: 
% Lmodel:   model struct of a linear problem
% solverset: at least contains field 'thesolver', which specify which
%            solver is used for solving fba problem, glpk as default. 

% Output:
% sol_lp: solution of the linear problem
% jingyi cai 2017-7

    switch solverset.thesolver
       case 'glpk'
            c=Lmodel.Model.obj;
            a=[Lmodel.Model.A;-Lmodel.Model.A];
            b=[Lmodel.Model.rhs;-Lmodel.Model.lhs];
            lb=Lmodel.Model.lb;
            ub=Lmodel.Model.ub;
            ctype=char('U'*ones(1,length(b)));
            vartype=char('C'*ones(1,length(c)));
            switch Lmodel.Model.sense
                case 'maximize'
                   thesense= -1;
                case 'minimize'
                   thesense =1;
            end
            [sol_lp.x, sol_lp.objval, sol_lp.status]=glpk (c, a, b, lb, ub, ctype, vartype,thesense);

            sol_lp.solver='glpk';
        case 'cplex'
            theLP=Cplex();
            theLP.Model=Lmodel.Model;
            sol_lp=theLP.solve();
            sol_lp.solver='cplex';

        case 'opti'
            switch Lmodel.Model.sense
                case 'maximize'
                   thesense= -1;
                case 'minimize'
                   thesense =1;
            end
          
            c=Lmodel.Model.obj;
            A=Lmodel.Model.A;
            % Setup Options
            opts = optiset('solver','clp'); %CLP solver
            rhs=Lmodel.Model.rhs;
            lhs=Lmodel.Model.lhs;
            lb=Lmodel.Model.lb;
            ub=Lmodel.Model.ub;
            theopt=opti('sense',thesense,'f',c,'lin',A,lhs,rhs,'bounds',lb,ub,'options',opts);
            [sol_lp.x, sol_lp.objval, sol_lp.status]=solve(theopt);
    end
end
function [merged_struct] = CombineStructs(struct_a,struct_b)

%%if one of the structres is empty do not combine
if isempty(struct_a)
    merged_struct=struct_b;
    return
end
if isempty(struct_b)
    merged_struct=struct_a;
    return
end
%%insert struct a
merged_struct=struct_a;
%%insert struct b
size_a=length(merged_struct);

f = fieldnames(struct_b);
for i = 1:length(f)
    merged_struct(size_a).(f{i}) = struct_b.(f{i});
end

end