function [sol, info] = exozymeLinear(X, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, findAllNE, params)
% Find a Nash equilibrium given the relative abundance of mutants 
% for the following simple exoenzyme game
%                      _________________
%                     |vU    vB         |
% Sucrose -> Glucose -|--> G --> biomass|
%        /|\          |   | \           |
%         |           |vE | _\|         |
%      exoenzyme <----|----   ATPM      |
%                     |_________________|    
%
% All fluxes (v) here are scaled by the biomass/relative abundance X
%
% d[exoenzyme]/dt = vE - degradationConstant * [exoenzyme] = 0
% [exoenzyme] = vE / degrad.Constant
% exoenzyme flux = glucose release = k_cat * [exoenzyme] * Sucrose / (Km + Sucrose) = k * vE
% Assume the cell can take all glucose available (therefore no the convexity observed in that paper)
% vU <= k * vE
% In a community setting, with the parameter epsCap for capture efficiency
% the glucose organism n is able to get = epsCap * k * v^n_E + X^n * (1 - epsCap) * k * sum(v^m_E for all m)
%
% For each member
% max v^n_B
% s.t.
%   v^n_U - v^n_B - costExo * v^n_E = ATPM * X^n (mass balance for intracellular glucose)
%   v^n_U <= epsCap * k * v^n_E + X^n * (1 - epsCap) * k * sum(v^m_E for all m)
%   v^n_U <= vUub * X^n
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
%     fixMutantVE a vector with a fixed value of exoenzyme production flux (vE) for each mutant. 
%                  Use NaN for mutants not to be fixed.
%     outObj      outer-level objective coefficients for the biomass
%                 production term (N x 1 vector for N mutants)
%     findAllNE   true to find all NE using integer cuts
%     params      cobra parameter structure
%
% OUTPUTS:
%    sol          cobra solution(s); structure (array)
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
    vUub = 1;  % maximum specific glucose uptake rate
end
if nargin < 6
    atpm = 0;  % quadratic cost for exoenzyme
end
if nargin < 7
    fixMutantVE = NaN(numel(X), 1);  % given exoenzyme production rate (vE) for mutants
elseif ~isempty(fixMutantVE) && ~(isscalar(fixMutantVE) && (~fixMutantVE || isnan(fixMutant)))
    if numel(fixMutantVE) ~= numel(X)
        error('fixMutantVE must be a vector with a fixed value of vE for each mutant. Use NaN for mutants not to be fixed')
    end
end
if nargin < 8
    outObj = 1 ./ X(:);  % quadratic cost for exoenzyme
end
if nargin < 9
    findAllNE = true;  % find all NE using integer Cuts
end

N = numel(X);
bigM = 1000; % big M

% uptake  biomass exoenzyme, dual variables, binary vars for complementary slackness
varList = {'vU', 'vB', 'vE', 'l', 'g', 'mu', 'ag', 'au', 'ab', 'ae', 'aeOn', 'amu'};
vartype = {'C',  'C',  'C',  'C', 'C', 'C',  'B',  'B',  'B',  'B',  'B',    'B'};
lb = [     0,    0,    0,-bigM*10, 0,   0,    0,   0,    0,    0,    0,      0];
ub = [10 * bigM * ones(1, 6), ones(1, 6)];

conList = {'mb', 'eflux', 'efluxCS', 'vUub', 'vUubCS', 'vUd', 'vUdCS', 'vBd', 'vBdCS', 'vEd', 'vEdCS', ...
    'gCS', 'vUCS', 'vBCS', 'vECS', 'muCS', 'vEon', 'vEoff'};

LP = struct();
LP.vartype = '';
LP.lb = zeros(0, 1);
LP.ub = zeros(0, 1);
v = struct();
c = struct();
colname = {};
rowname = {};
for j = 1:numel(varList)
    v.(varList{j}) = (N * (j - 1) + 1):(N * j);
    LP.vartype(v.(varList{j})) = vartype{j};
    LP.lb(v.(varList{j})) = lb(j);
    LP.ub(v.(varList{j})) = ub(j);
    colname(v.(varList{j})) = strcat(varList{j}, '_', strtrim(cellstr(num2str((1:N)'))));
end
for j = 1:numel(conList)
    c.(conList{j}) = (N * (j - 1) + 1):(N * j);
    rowname(c.(conList{j})) = strcat(conList{j}, '_', strtrim(cellstr(num2str((1:N)'))));
end
nV = v.(varList{end})(end);
nC = c.(conList{end})(end);

[row, col, e] = deal([]);
LP.b = zeros(0, 1);
LP.csense = '';
% mass balance for intracellular sugar
% v_u - v_b - costExo * v_E = atpm * X
row = [row, c.mb, c.mb, c.mb];
col = [col, v.vU, v.vB, v.vE];
e = [e, ones(1, N), -ones(1, N), -costExo * ones(1, N)];
LP.b(c.mb) = atpm * X(:);
LP.csense(c.mb) = 'E';

% v^n_u - k*(epsCap + X^n*(1 - epsCap)) * v^n_E - X^n*k*(1-epsCap)*sum(m!=n, v^m_E) <= 0
row = [row, repmat(c.eflux, 1, N + 1)];
col = [col, v.vU];
for j = 1:N
    col = [col, v.vE([j:N, 1:(j - 1)])];
end
e = [e, ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), repmat(-X(:)' * k * (1 - epsCap), 1, N - 1)];
LP.b(c.eflux) = 0;
LP.csense(c.eflux) = 'L';

% v^n_u - k(eps+X(1-eps)*v^n_E - X^n*k*(1-epsCap)*sum(m!=n, v^m_E) - M*a_g >= -M
row = [row, repmat(c.efluxCS, 1, N + 2)];
col = [col, v.vU];
for j = 1:N
    col = [col, v.vE([j:N, 1:(j - 1)])];
end
col = [col, v.ag];
e = [e, ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), repmat(-X(:)' * k * (1 - epsCap), 1, N - 1), -bigM * ones(1, N)];
LP.b(c.efluxCS) = -bigM;
LP.csense(c.efluxCS) = 'G';

% v^n_u <= vUub * X
row = [row, c.vUub];
col = [col, v.vU];
e = [e, ones(1, N)];
LP.b(c.vUub) = vUub * X;
LP.csense(c.vUub) = 'L';

% v^n_u - M*amu >= vUub * X - M
row = [row, c.vUubCS, c.vUubCS];
col = [col, v.vU, v.amu];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vUubCS) = vUub * X - bigM;
LP.csense(c.vUubCS) = 'G';

% l + g + mu >= 0
row = [row, repmat(c.vUd, 1, 3)];
col = [col, v.l, v.g, v.mu];
e = [e, ones(1, N * 3)];
LP.b(c.vUd) = 0;
LP.csense(c.vUd) = 'G';

% l + g + mu + M*a_u <= M
row = [row, repmat(c.vUdCS, 1, 4)];
col = [col, v.l, v.g, v.mu, v.au];
e = [e, ones(1, N * 3), bigM * ones(1, N)];
LP.b(c.vUdCS) = bigM;
LP.csense(c.vUdCS) = 'L';

% l <= -1
row = [row, c.vBd];
col = [col, v.l];
e = [e, ones(1, N)];
LP.b(c.vBd) = -1;
LP.csense(c.vBd) = 'L';

% l - M*a_b >= -1 - M
row = [row, c.vBdCS, c.vBdCS];
col = [col, v.l, v.ab];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vBdCS) = -1 - bigM;
LP.csense(c.vBdCS) = 'G';

% -costExo * l - k(eps+X(1-eps)*g >= 0
row = [row, c.vEd, c.vEd];
col = [col, v.l, v.g];
e = [e, -costExo * ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap))];
LP.b(c.vEd) = 0;
LP.csense(c.vEd) = 'G';

% -costExo*l - k(eps+X(1-eps)*g + M*a_E <= M
row = [row, repmat(c.vEdCS, 1, 3)];
col = [col, v.l, v.g, v.ae];
e = [e, -costExo * ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), bigM * ones(1, N)];
LP.b(c.vEdCS) = bigM;
LP.csense(c.vEdCS) = 'L';

% g - M*a_g <= 0
row = [row, repmat(c.gCS, 1, 2)];
col = [col, v.g, v.ag];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.gCS) = 0;
LP.csense(c.gCS) = 'L';

% vU - M*a_u <= 0
row = [row, repmat(c.vUCS, 1, 2)];
col = [col, v.vU, v.au];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vUCS) = 0;
LP.csense(c.vUCS) = 'L';

% vB - M*a_b <= 0
row = [row, repmat(c.vBCS, 1, 2)];
col = [col, v.vB, v.ab];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vBCS) = 0;
LP.csense(c.vBCS) = 'L';

% vE - M*a_i <= 0
row = [row, repmat(c.vECS, 1, 2)];
col = [col, v.vE, v.ae];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vECS) = 0;
LP.csense(c.vECS) = 'L';

% mu - M*amu <= 0
row = [row, repmat(c.muCS, 1, 2)];
col = [col, v.mu, v.amu];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.muCS) = 0;
LP.csense(c.muCS) = 'L';

% vE - M*aeOn >= -M + epsFlux
epsFlux = 1e-3;
row = [row, repmat(c.vEon, 1, 2)];
col = [col, v.vE, v.aeOn];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vEon) = -bigM + epsFlux;
LP.csense(c.vEon) = 'G';

% vE - M*aeOn <= 0
row = [row, repmat(c.vEoff, 1, 2)];
col = [col, v.vE, v.aeOn];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vEoff) = 0;
LP.csense(c.vEoff) = 'L';

LP.A = sparse(row, col, e, nC, nV);
LP.c = zeros(nV, 1);
LP.c(v.vB) = outObj;
LP.osense = -1;
LP.x0 = [];

if ~(isscalar(fixMutantVE) && (~fixMutantVE || isnan(fixMutant)))
    mutantIDs = ~isnan(fixMutantVE);
    [LP.lb(v.vE(mutantIDs)), LP.ub(v.vE(mutantIDs))] = deal(fixMutantVE(mutantIDs));
    %     % relaxing the two constraints below will find solutions that are not NE
    %     LP.rhs(c.vEd(mutantIDs)) = -inf;
    %     LP.rhs(c.vEdCS(mutantIDs)) = inf;
end
if exist('params', 'var')
    sol = solveCobraMILP(LP, params);
else
    sol = solveCobraMILP(LP);
end
if findAllNE
    c.intCut = [];
    nIC = 0;
    solCur = sol;
    while solCur.stat == 1
        nIC = nIC + 1;
        sol(nIC) = solCur;
        % integer cut
        % sum(aeOn_prev == 1, aeOn) - sum(aeOn_prev == 0, aeOn) <= sum(aeOn_prev) - 1
        LP.A = [LP.A; ...
            sparse(1, v.aeOn, (solCur.full(v.aeOn) > 0.9) - (solCur.full(v.aeOn) < 0.1), 1, nV)];
        LP.b(end + 1) = sum(solCur.full(v.aeOn) > 0.9) - 1;
        LP.csense(end + 1) = 'L';
        c.intCut(end + 1) = size(LP.A, 1);
        rowname(end + 1) = {['integerCut' num2str(nIC)]};
        if exist('params', 'var')
            solCur = solveCobraMILP(LP, params);
        else
            solCur = solveCobraMILP(LP);
        end
    end
end

info = struct();
info.v = v;
info.c = c;
info.LP = LP;
info.colname = colname;
info.rowname = rowname;
