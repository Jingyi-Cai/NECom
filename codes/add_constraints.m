function Lmodel=add_constraints(Lmodel,lhs,A,rhs,connames)
%{
Function: Add constraints to linear models: lhs<=AX<=rhs
Input: 
    Lmodel: model struct with at least 'Model' field.
    lhs: M by 1 real number vector, left-hand-side of constraints 
    rhs: M by 1 real number vector, right-hand-side of constraints 
    A: M by N matrix, the linear model matrix
    connames: M by 1 string cells,  constraints names( optional)     
Output: 
    Lmodel: model struct with at least 'Model' field.
Jingyi Cai 2017-07-31
%}



if nargin<5
    connames=[];
end

len1=size(Lmodel.Model.A,2);
len2=size(Lmodel.Model.lb,1);
Lmodel.Model.A=[Lmodel.Model.A,zeros(size(Lmodel.Model.A,1),len2-len1)];
Lmodel.Model.A=[Lmodel.Model.A;A];
Lmodel.Model.lhs=[Lmodel.Model.lhs;lhs];
Lmodel.Model.rhs=[Lmodel.Model.rhs;rhs];


if isempty(connames)
    connames=char(ones(length(lhs),1)*'na');
end

Lmodel.Model.conname=[Lmodel.Model.conname;cellstr(connames)];

end