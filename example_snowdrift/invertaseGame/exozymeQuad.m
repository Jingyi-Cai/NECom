function [sol, info] = exozymeQuad(X, epsCap, costExo, k, r, vUub, qc, fixMutantVE, outObj, findAllNE, params)
% Find a Nash equilibrium given the relative abundance of mutants 
% for the following simple exoenzyme game
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
%   v^n_ATPM >= (r / v^max) * (v^n_B + qc*v^n_E) ^ 2  (ATPM changing quadratically with biomass and exoenzyme production)
%   v^n_U, v^n_B, v^n_E >= 0
%
% INPUTS:
%     X           relative abundance profile (N x 1 vector for N mutants)
%     epsCap      capture efficiency, percentage of sugar that exoenzyme producer 
%                 will have advantage for uptake 0 <= epsCap <= 1
%     costExo     sugar per exoenzyme production needed (equivalent to cost)
%     k           k_effective for exoenzyme divided by the degradation rate constant
%                 assuming steady state for the investase in the media
%     r           scaling factor for ATPM dependence on biomass and exoenzyme production
%     vUub        maximum sugar uptake rate
%     qc          scaling factor associated with the quadratic terms
%                 in the ATPM requirement involving exoenzyme production
%     fixMutantVE a vector with a fixed value of exoenzyme production flux (vE) for each mutant. 
%                  Use NaN for mutants not to be fixed.
%     outObj      outer-level objective coefficients for the biomass
%                 production term (N x 1 vector for N mutants)
%     findAllNE   true to find all NE using integer cuts
%     params      gurobi parameter structure
%
% OUTPUTS:
%    sol          gurobi solution structure
%    info         structure with the following:
%                 *.LP       cobra LP problem structure
%                 *.c        indices for linear constraints
%                 *.v        indices for variables
%                 *.rowname  cell array of constraint names
%                 *.colname  cell array of variable names

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
    r = 0.1;  % scaling factor for GAM dependence on biomass production
end
if nargin < 6
    vUub = 1;  % maximum specific glucose uptake rate
end
if nargin < 7
    qc = 1;  % quadratic cost for exoenzyme
end
if nargin < 8
    fixMutantVE = NaN(numel(X), 1);  % given exoenzyme production rate (vE) for mutants
elseif ~isempty(fixMutantVE) && ~(isscalar(fixMutantVE) && (~fixMutantVE || isnan(fixMutant)))
    if numel(fixMutantVE) ~= numel(X)
        error('fixMutantVE must be a vector with a fixed value of vE for each mutant. Use NaN for mutants not to be fixed')
    end
end
if nargin < 9
    outObj = 1 ./ X(:);  % quadratic cost for exoenzyme
end
if nargin < 10
    findAllNE = true;  % find all NE using integer Cuts
end

N = numel(X);
bigM = 1000; % big M

objCoeff = ones(N,1);%(1 - X) ./ X;

%       uptake  biomass exoenzyme, dual variables, binary vars for comp. slack
varList = {'vU', 'vB', 'vE', 'vM', 'l',   'g', 'mu', 'alpha', 'ag', 'au', 'ae', 'amu', 'ab', 'aeOn'};
vartype = {'C',  'C',  'C',  'C',  'C',   'C', 'C',  'C',     'B',  'B',  'B',  'B',   'B',  'B'};
lb = [     0,     0,    0,    0, -bigM*10, 0,   0,    0,       0,    0,    0,    0,    0,    0];
ub = [10 * bigM * ones(1, 8), ones(1, 6)];

linconList = {'mb', 'eflux', 'efluxCS', 'vUub', 'vUubCS', 'vUd', 'vUdCS', 'vEd', 'vEdCS', ...
    'vMd', 'gCS', 'muCS', 'vUCS', 'vECS', 'vBCS', 'vEon', 'vEoff'};

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
% mass balance for intracellular glucose
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

% v^n_u - k(eps+X(1-eps)*v^n_E - X^n*k*(1-epsCap)*sum(m!=n, v^m_E) - M*a_g >= -M
row = [row, repmat(c.efluxCS, 1, N + 2)];
col = [col, v.vU];
for j = 1:N
    col = [col, v.vE([j:N, 1:(j - 1)])];
end
col = [col, v.ag];
e = [e, ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), repmat(-X(:)' * k * (1 - epsCap), 1, N - 1), -bigM * ones(1, N)];
LP.rhs(c.efluxCS) = -bigM;
LP.sense(c.efluxCS) = '>';

% v^n_u <= vUub * X
row = [row, c.vUub];
col = [col, v.vU];
e = [e, ones(1, N)];
LP.rhs(c.vUub) = vUub * X;
LP.sense(c.vUub) = '<';

% v^n_u - M*amu >= vUub * X - M
row = [row, c.vUubCS, c.vUubCS];
col = [col, v.vU, v.amu];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.vUubCS) = vUub * X - bigM;
LP.sense(c.vUubCS) = '>';

% l + g + mu >= 0
row = [row, repmat(c.vUd, 1, 3)];
col = [col, v.l, v.g, v.mu];
e = [e, ones(1, N * 3)];
LP.rhs(c.vUd) = 0;
LP.sense(c.vUd) = '>';

% l + g + mu + M*a_u <= M
row = [row, repmat(c.vUdCS, 1, 4)];
col = [col, v.l, v.g, v.mu, v.au];
e = [e, ones(1, N * 3), bigM * ones(1, N)];
LP.rhs(c.vUdCS) = bigM;
LP.sense(c.vUdCS) = '<';

% % l = -1
% row = [row, c.vBd];
% col = [col, v.l];
% e = [e, ones(1, N)];
% LP.rhs(c.vBd) = -1;
% LP.sense(c.vBd) = 'E';

% % l - M*a_b >= -1 - M
% row = [row, c.vBdCS, c.vBdCS];
% col = [col, v.l, v.ab];
% e = [e, ones(1, N), -bigM * ones(1, N)];
% LP.rhs(c.vBdCS) = -1 - bigM;
% LP.sense(c.vBdCS) = 'G';

% (qc - costExo) * l - k(eps+X(1-eps)*g >= -qc
row = [row, c.vEd, c.vEd];
col = [col, v.l, v.g];
e = [e, (qc - costExo) * ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap))];
LP.rhs(c.vEd) = -qc * objCoeff;
LP.sense(c.vEd) = '>';

% (qc-costExo)*l - k(eps+X(1-eps)*g + M*a_E <= M - qc * objCoeff
row = [row, repmat(c.vEdCS, 1, 3)];
col = [col, v.l, v.g, v.ae];
e = [e, (qc - costExo) * ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), bigM * ones(1, N)];
LP.rhs(c.vEdCS) = bigM - qc * objCoeff;
LP.sense(c.vEdCS) = '<';

% g - M*a_g <= 0
row = [row, repmat(c.gCS, 1, 2)];
col = [col, v.g, v.ag];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.gCS) = 0;
LP.sense(c.gCS) = '<';

% mu - M*amu <= 0
row = [row, repmat(c.muCS, 1, 2)];
col = [col, v.mu, v.amu];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.muCS) = 0;
LP.sense(c.muCS) = '<';

% vU - M*a_u <= 0
row = [row, repmat(c.vUCS, 1, 2)];
col = [col, v.vU, v.au];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.vUCS) = 0;
LP.sense(c.vUCS) = '<';

% % vB - M*a_b <= 0
% row = [row, repmat(c.vBCS, 1, 2)];
% col = [col, v.vB, v.ab];
% e = [e, ones(1, N), -bigM * ones(1, N)];
% LP.rhs(c.vBCS) = 0;
% LP.sense(c.vBCS) = 'L';

% vE - M*a_i <= 0
row = [row, repmat(c.vECS, 1, 2)];
col = [col, v.vE, v.ae];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.vECS) = 0;
LP.sense(c.vECS) = '<';

% -l - alpha = 0
row = [row, repmat(c.vMd, 1, 2)];
col = [col, v.l, v.alpha];
e = [e, -ones(1, N * 2)];
LP.rhs(c.vMd) = 0;
LP.sense(c.vMd) = '=';

% vB - M*a_b <= 0
row = [row, repmat(c.vBCS, 1, 2)];
col = [col, v.vB, v.ab];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.vBCS) = 0;
LP.sense(c.vBCS) = '<';

% vE - M*aeOn >= -M + epsFlux
epsFlux = 1e-3;
row = [row, repmat(c.vEon, 1, 2)];
col = [col, v.vE, v.aeOn];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.vEon) = -bigM + epsFlux;
LP.sense(c.vEon) = '>';

% vE - M*aeOn <= 0
row = [row, repmat(c.vEoff, 1, 2)];
col = [col, v.vE, v.aeOn];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.rhs(c.vEoff) = 0;
LP.sense(c.vEoff) = '<';

% (r/vU^max/X) * (vB + qc*vE)^2 - vM = 0 (alpha > 0 and lambda < 0, otherwise infeasible)
nQC = 0;
for j = 1:N
    nQC = nQC + 1;
    LP.quadcon(nQC).Qrow = [v.vB(j); v.vB(j); v.vE(j); v.vE(j)];
    LP.quadcon(nQC).Qcol = [v.vB(j); v.vE(j); v.vB(j); v.vE(j)];
    LP.quadcon(nQC).Qval = (r / vUub / X(j)) * [1; qc; qc; qc ^ 2];
%     LP.quadcon(2 * j - 1).Qrow = [v.vB(j); v.vE(j)];
%     LP.quadcon(2 * j - 1).Qcol = [v.vB(j); v.vE(j)];
    %LP.quadcon(2 * j - 1).Qval = (r / vUub / X(j))* ones(4, 1);

    LP.quadcon(nQC).q = full(sparse(v.vM(j), 1, -1, nV,1));
    LP.quadcon(nQC).rhs = 0;
    LP.quadcon(nQC).sense = '=';
    LP.quadcon(nQC).name = ['vMlb_' num2str(j)];
    
    % (2r/vU/X) * alpha * (vB + qc*vE) - l >= 1 
    nQC = nQC + 1;
    LP.quadcon(nQC).Qrow = [v.alpha(j); v.vB(j); v.alpha(j); v.vE(j)];
    LP.quadcon(nQC).Qcol = [v.vB(j); v.alpha(j); v.vE(j); v.alpha(j)];
    LP.quadcon(nQC).Qval = (r / vUub / X(j)) * qc * [1; 1; qc; qc];
    LP.quadcon(nQC).q = full(sparse(v.l(j), 1, -1, nV,1));
    LP.quadcon(nQC).rhs = objCoeff(j);
    LP.quadcon(nQC).sense = '>';
    LP.quadcon(nQC).name = ['vBd_' num2str(j)];
    
    % (2r/vU/X) * alpha * (vB + qc*vE) - l + bigM * a_b <= bigM + 1 
    nQC = nQC + 1;
    LP.quadcon(nQC).Qrow = [v.alpha(j); v.vB(j); v.alpha(j); v.vE(j)];
    LP.quadcon(nQC).Qcol = [v.vB(j); v.alpha(j); v.vE(j); v.alpha(j)];
    LP.quadcon(nQC).Qval = (r / vUub / X(j)) * qc * [1; 1; qc; qc];
    LP.quadcon(nQC).q = full(sparse([v.l(j); v.ab(j)], 1, [-1; bigM], nV,1));
    LP.quadcon(nQC).rhs = bigM + objCoeff(j);
    LP.quadcon(nQC).sense = '<';
    LP.quadcon(nQC).name = ['vBdCS_' num2str(j)];
end


LP.A = sparse(row, col, e, nC, nV);
LP.obj = zeros(nV, 1);
LP.obj(v.vB) = outObj;
LP.modelsense  = 'max';
LP.varnames = colname;
LP.constrnames = rowname;

if ~(isscalar(fixMutantVE) && (~fixMutantVE || isnan(fixMutant)))
    mutantIDs = ~isnan(fixMutantVE);
    [LP.lb(v.vE(mutantIDs)), LP.ub(v.vE(mutantIDs))] = deal(fixMutantVE(mutantIDs));
    %     % relaxing the two constraints below will find solutions that are not NE
    %     LP.rhs(c.vEd(mutantIDs)) = -inf;
    %     LP.rhs(c.vEdCS(mutantIDs)) = inf;
end


params.NonConvex = 2;
sol = gurobi(LP, params);

if findAllNE
    c.intCut = [];
    nIC = 0;
    solCur = sol;
    while strcmp(solCur.status, 'OPTIMAL')
        nIC = nIC + 1;
        sol(nIC) = solCur;
        % integer cut
        % sum(aeOn_prev == 1, aeOn) - sum(aeOn_prev == 0, aeOn) <= sum(aeOn_prev) - 1
        LP.A = [LP.A; ...
            sparse(1, v.aeOn, (solCur.x(v.aeOn) > 0.9) - (solCur.x(v.aeOn) < 0.1), 1, nV)];
        LP.rhs(end + 1) = sum(solCur.x(v.aeOn) > 0.9) - 1;
        LP.sense(end + 1) = '<';
        LP.constrnames(end + 1) = {['integerCut' num2str(nIC)]};
        c.intCut(end + 1) = size(LP.A, 1);
        rowname(end + 1) = {['integerCut' num2str(nIC)]};
        solCur = gurobi(LP, params);
    end
end

info = struct();
info.v = v;
info.c = c;
info.LP = LP;
info.colname = colname;
info.rowname = rowname;
