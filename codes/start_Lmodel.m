function Lmodel=start_Lmodel(solverset)
% FUNCTION: start building a linear model.
% Jingyi Cai, July 2016
        Lmodel.Model.A=[];
        Lmodel.Model.lhs=[];
        Lmodel.Model.rhs=[];
        Lmodel.Model.obj=[];
        Lmodel.Model.lb=[];
        Lmodel.Model.ub=[];
        Lmodel.Model.vtype=[];  
        Lmodel.Model.conname='';
        Lmodel.Model.varname='';
if exist('solverset')
    
    Lmodel =solver_setup_proper(Lmodel,solverset.solverValues); 

end               

end


function lp=solver_setup_proper(lp,values)
% solver_setup_proper setup CPLEX tolerance parameters properly
% Jingyi Cai, July 2016
%lp=Cplex();
% --------------tolerances-------------------
%feasibility tolerance http://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.2/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpRHS.html
%the degree to which the basic variables of a model may violate their bounds
%min 1e-9, 1e-6 default
lp.Param.simplex.tolerances.feasibility.Cur=values(1);

%convergence tolerance http://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.2/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarEpComp.html
%Sets the tolerance on complementarity for convergence.
%min 1e-12, 1e-8 default
lp.Param.barrier.convergetol.Cur=values(2);
%integrality tolerance http://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.2/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpInt.html
%Specifies the amount by which an integer variable can be different from an integer and still be considered feasible.
%min 0, 1e-5 default
lp.Param.mip.tolerances.integrality.Cur=values(3);

%relative/absolute MIP gap tolerance http://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.2/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpAGap.html
%Sets an relative/absolute tolerance on the gap between the best integer objective and the objective of the best node remaining.
%min 0 , 1e-6 default 
%min 0 , 1e-4 default
lp.Param.mip.tolerances.absmipgap.Cur=1e-8;
lp.Param.mip.tolerances.mipgap.Cur=1e-8;
%optimality tolerance http://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.2/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpOpt.html
%Influences the reduced-cost tolerance for optimality.
%min 1e-9, le-6 default
lp.Param.simplex.tolerances.optimality.Cur=1e-9;

end

