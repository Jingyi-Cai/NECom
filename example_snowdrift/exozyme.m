% A simple exoenzyme game
%                      _________________
%                     |vU    vB         |
% Sucrose -> Glucose -|--> G --> biomass|
%        /|\          |   | \           |
%         |           |vI | _\|         |
%      Invertase <----|----   ATPM      |
%                     |_________________|    
%
% All fluxes (v) here are scaled by the biomass/relative abundance X
%
% d[Invertase]/dt = vI - degradationConstant * [Invertase] = 0
% [Invertase] = vI / degrad.Constant
% Invertase flux = glucose release = k_cat * [Invertase] * Sucrose / (Km + Sucrose) = k * vI
% Assume the cell can take all glucose available (therefore no the convexity observed in that paper)
% vU <= k * vI
% In a community setting, with the parameter epsCap for capture efficiency
% the glucose organism n is able to get = epsCap * k * v^n_I + X^n * (1 - epsCap) * k * sum(v^m_I for all m)
%
% For each member
% max v^n_B
% s.t.
%   v^n_U - v^n_B - v^n_I = ATPM * X^n (mass balance for intracellular glucose)
%   v^n_U <= epsCap * k * v^n_I + X^n * (1 - epsCap) * k * sum(v^m_I for all m)
%   v^n_U <= vUub * X^n
%   v^n_U, v^n_B, v^n_I >= 0

N = 10;  % N: number of mutants

% epsCap: capture efficiency 
% (advantage for invertase producer, 0 = none, 1 = get all it breaks down)
epsCap = 0.0;

% A pseudo-kinetic parameter
% k_effective for invertase divided by the degradation rate constant
% assuming steady state for the investase in the media
k = 10;

X = repmat(1/N, N, 1);  % relative abundance
atpm = 0.1;  % atpm: NGAM value
bigM = 1000; % big M
vUub = 1;  % maximum specific glucose uptake rate
%       uptake  biomass invertase, dual variables, binary vars for comp. slack
varList = {'vU', 'vB', 'vI', 'l', 'g', 'mu', 'ag', 'au', 'ab', 'ai', 'amu'};
vartype = {'C',  'C',  'C',  'C', 'C', 'C',  'B',  'B',  'B',  'B',  'B'};
lb = [     0,    0,    0,-bigM*10, 0,   0,    0,   0,    0,    0,    0];
ub = [10 * bigM * ones(1, 6), ones(1, 5)];

conList = {'mb', 'eflux', 'efluxCS', 'vUub', 'vUubCS', 'vUd', 'vUdCS', 'vBd', 'vBdCS', 'vId', 'vIdCS', ...
    'gCS', 'vUCS', 'vBCS', 'vICS', 'muCS'};

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
% mass balance for intracellular glucose
% v_u - v_b - v_I = atpm * X
row = [row, c.mb, c.mb, c.mb];
col = [col, v.vU, v.vB, v.vI];
e = [e, ones(1, N), -ones(1, N), -ones(1, N)];
LP.b(c.mb) = atpm * X(:);
LP.csense(c.mb) = 'E';

% v^n_u - k*(epsCap + X^n*(1 - epsCap)) * v^n_I - X^n*k*(1-epsCap)*sum(m!=n, v^m_I) <= 0
row = [row, repmat(c.eflux, 1, N + 1)];
col = [col, v.vU];
for j = 1:N
    col = [col, v.vI([j:N, 1:(j - 1)])];
end
e = [e, ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), repmat(-X(:)' * k * (1 - epsCap), 1, N - 1)];
LP.b(c.eflux) = 0;
LP.csense(c.eflux) = 'L';

% v^n_u - k(eps+X(1-eps)*v^n_I - X^n*k*(1-epsCap)*sum(m!=n, v^m_I) - M*a_g >= -M
row = [row, repmat(c.efluxCS, 1, N + 2)];
col = [col, v.vU];
for j = 1:N
    col = [col, v.vI([j:N, 1:(j - 1)])];
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

% l + g + mu>= 0
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

% -l - k(eps+X(1-eps)*g >= 0
row = [row, c.vId, c.vId];
col = [col, v.l, v.g];
e = [e, -ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap))];
LP.b(c.vId) = 0;
LP.csense(c.vId) = 'G';

% -l - k(eps+X(1-eps)*g + M*a_I <= M
row = [row, repmat(c.vIdCS, 1, 3)];
col = [col, v.l, v.g, v.ai];
e = [e, -ones(1, N), -k * (epsCap + X(:)' * (1 - epsCap)), bigM * ones(1, N)];
LP.b(c.vIdCS) = bigM;
LP.csense(c.vIdCS) = 'L';

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

% vI - M*a_i <= 0
row = [row, repmat(c.vICS, 1, 2)];
col = [col, v.vI, v.ai];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.vICS) = 0;
LP.csense(c.vICS) = 'L';

% mu - M*amu <= 0
row = [row, repmat(c.muCS, 1, 2)];
col = [col, v.mu, v.amu];
e = [e, ones(1, N), -bigM * ones(1, N)];
LP.b(c.muCS) = 0;
LP.csense(c.muCS) = 'L';

LP.A = sparse(row, col, e, nC, nV);
LP.c = zeros(nV, 1);
LP.c(v.vB) = 1;
LP.osense = -1;
LP.x0 = [];


sol = solveCobraMILP(LP);
if sol.stat == 1
    printModel(LP, [], 'obj', 'col', colname, 'row', rowname,'sol',sol)
    printModel(LP, [c.mb,c.eflux], 'con', 'col', colname, 'row', rowname,'sol',sol)
else
    fprintf('Seems no optimal solution found!\n')
end