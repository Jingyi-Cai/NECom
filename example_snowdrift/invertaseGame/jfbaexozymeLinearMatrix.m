function [gr, info] = jfbaexozymeLinearMatrix(X, epsCap, costExo, k, atpm, vUub, params)
% Compute a payoff matrix for the following simple exoenzyme game
%                      _________________
%                     |vU    vB         |
% Sucrose -> Glucose -|--> G --> biomass|
%        /|\          |   | \           |
%         |           |vE | _\|         |
%      Exoenzyme <----|----   ATPM      |
%                     |_________________|    
%
% All fluxes (v) here are scaled by the biomass, i.e., already multiplied by the relative abundance X
%
% d[Exoenzyme]/dt = vE - degradationConstant * [Exoenzyme] = 0
% [Exoenzyme] = vE / degradationConstant
% Exoenzyme flux = glucose release = k_cat * [Exoenzyme] * Sucrose / (Km + Sucrose) = k * vE
% Assume the cell can take all sugar available vU <= k * vE
% In a community setting, with the parameter epsCap for capture efficiency
% of the sugar, organism n is able to get:
%     epsCap * k * v^n_E + X^n * (1 - epsCap) * k * sum(v^m_E for all m)
% Additionally, the maintenance requirement is modified to change to
% increase with the biomass production to make it a convex function of
% uptake
% 
% For each member
% max v^n_B
% s.t.
%   v^n_U - v^n_B - v^n_E - v^n_M * X^n = 0 (mass balance for intracellular sugar)
%   v^n_U <= epsCap * k * v^n_E + X^n * (1 - epsCap) * k * sum(v^m_E for all m) (sugar available)
%   v^n_U <= vUub * X^n (max. uptake capacity)
%   v^n_M >= atpm  (ATPM changing quadratically with biomass and exoenzyme production)
%   v^n_U, v^n_B, v^n_E >= 0
%
% INPUTS:
%     X           relative abundance profile (N x 1 vector for N mutants)
%     epsCap      capture efficiency, percentage of sugar that exoenzyme producer 
%                 will have advantage for uptake 0 <= epsCap <= 1
%     costExo     sugar per exoenzyme production needed (equivalent to cost)
%     k           k_effective for exoenzyme divided by the degradation rate constant
%                 assuming steady state for the investase in the media
%     atpm        constant ATPM requirement (NGAM)
%     vUub        maximum sugar uptake rate
%     nStep       number of steps to slice vE, or a vector of values to
%                 test directly
%     params      gurobi parameter structure
%
% OUTPUTS:
%    gr           [nStep^N by N matrix] growth rates in each row are the growth rates for all mutants
%                 in each of the nStep^N simulated tuples of (vE_1, vE_2, ..., vE_N)
%    vErange      [nStep by 1 vector] vE values simulated to generate growth rates
%    combinations [nStep^N by N matrix] vErange(combinations(k, :)) is the
%                 vE tuple simulated in the k-th simulation, 
%                 corresponding to the growth rates in gr(k, :)
%    info         structure with the following:
%                 *.LP       cobra LP problem structure
%                 *.c        indices for linear constraints
%                 *.v        indices for variables
%                 *.rowname  cell array of constraint names
%                 *.colname  cell array of variable names
%
% If N = 2, to generate a payoff matrix for mutant 1, for example, 
% use full(sparse(combinations(:, 1), combinations(:, 2), gr(:, 1))

if nargin < 1
    X = [0.5; 0.5];
end
if nargin < 2
    % epsCap: capture efficiency
    % (advantage for exoenzyme producer, 0 = none, 1 = get all it breaks down)
    epsCap = 0.01;
end
if nargin < 3
    % cost of exoenzyme
    costExo = 0.1;
end
if nargin < 4
    % A pseudo-kinetic parameter
    % k_effective for exoenzyme divided by the degradation rate constant
    % assuming steady state for the investase in the media
    k = 10;
end
if nargin < 5
    atpm = 0;  % scaling factor for GAM dependence on biomass production
end
if nargin < 6
    vUub = 1;  % maximum specific sugar uptake rate
end


N = numel(X);
bigM = 1000; % big M

%objCoeff = (1 - X) ./ X;

%       uptake  biomass exoenzyme, dual variables, binary vars for comp. slack
varList = {'vU', 'vB', 'vE', 'vM'};
vartype = {'C',  'C',  'C',  'C'};
lb = [     0, 0,  0,   0];
ub = 10 * bigM * ones(1, 4);

linconList = {'mb', 'eflux', 'vUub', 'vMlb'};

LP = struct();
LP.vtype = '';
LP.lb = zeros(0, 1);
LP.ub = zeros(0, 1);
v = struct();
c = struct();
colname = {};
rowname = {};
for j = 1:numel(varList)
    v.(varList{j}) = (N * (j - 1) + 1):(N * j);
    LP.vtype(v.(varList{j})) = vartype{j};
    LP.lb(v.(varList{j})) = lb(j);
    LP.ub(v.(varList{j})) = ub(j);
    colname(v.(varList{j})) = strcat(varList{j}, '_', strtrim(cellstr(num2str((1:N)'))));
end
for j = 1:numel(linconList)
    c.(linconList{j}) = (N * (j - 1) + 1):(N * j);
    rowname(c.(linconList{j})) = strcat(linconList{j}, '_', strtrim(cellstr(num2str((1:N)'))));
end
nV = v.(varList{end})(end);
nC = c.(linconList{end})(end);

[row, col, e] = deal([]);
LP.rhs = zeros(0, 1);
LP.sense = '';
% mass balance for intracellular sugar
% v_u - v_b - costExo * v_E - v_M = 0
row = [row, c.mb, c.mb, c.mb, c.mb];
col = [col, v.vU, v.vB, v.vE, v.vM];
e = [e, ones(1, N), -ones(1, N), -costExo * ones(1, N), -ones(1, N)];
LP.rhs(c.mb) = 0;
LP.sense(c.mb) = '=';

% v^n_u - k*(epsCap + X^n*(1 - epsCap)) * v^n_E - X^n*k*(1-epsCap)*sum(m!=n, v^m_E) <= 0
row = [row, repmat(c.eflux, 1, N + 1)];
col = [col, v.vU];
for j = 1:N
    col = [col, v.vE([j:N, 1:(j - 1)])];
end
e = [e, ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), repmat(-X(:)' * k * (1 - epsCap), 1, N - 1)];
LP.rhs(c.eflux) = 0;
LP.sense(c.eflux) = '<';

% v^n_u <= vUub * X
row = [row, c.vUub];
col = [col, v.vU];
e = [e, ones(1, N)];
LP.rhs(c.vUub) = vUub * X;
LP.sense(c.vUub) = '<';

% v^n_M >= atpm * X
row = [row, c.vMlb];
col = [col, v.vM];
e = [e, ones(1, N)];
LP.rhs(c.vMlb) = atpm * X;
LP.sense(c.vMlb) = '>';

LP.A = sparse(row, col, e, nC, nV);
LP.obj = zeros(nV, 1);
LP.varnames = colname;
LP.constrnames = rowname;
LP.modelsense  = 'max';



LP.obj(:) = 0;
LP.obj(v.vB) = 1./X;


sol = gurobi(LP, params);
if strcmp(sol.status, 'OPTIMAL')
    gr(1, :) = sol.x(v.vB)' ./ X(:)';
else
    gr(1, :) = NaN;
end



info = struct();
info.v = v;
info.c = c;
info.LP = LP;
info.colname = colname;
info.rowname = rowname;
