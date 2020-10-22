function Lmodel=add_variables(Lmodel,objvar,lb,ub,vtype,varnames)
%{
Function: Add variables to linear models
Input: 
    Lmodel: model struct with at least 'Model' field.
    objvar: N by 1 objective indicator vector, in which 0 for non-objective 1
    for objective
    lb: N by 1 real number vector, lower bounds of variable 
    ub: N by 1 real number vector, lower bounds of variable 
    vtype: 1 by N character string, 'C','I','B' variable type
    varnames: N by 1 string cells,  variable names( optional)     
Output: 
    Lmodel: model struct with at least 'Model' field.
Jingyi Cai 2017-07-31
%}

if nargin<5
    vtype=char(ones(1,length(lb))*'C');
end

if nargin<6
    varnames='';
end



Lmodel.Model.obj=[Lmodel.Model.obj;objvar];
Lmodel.Model.lb=[Lmodel.Model.lb;lb];
Lmodel.Model.ub=[Lmodel.Model.ub;ub];
Lmodel.Model.vtype=[Lmodel.Model.vtype,vtype];  

if isempty(varnames)
    varnames=char(ones(length(lb),1)*'na');
end


Lmodel.Model.varname=[Lmodel.Model.varname;cellstr(varnames)];
