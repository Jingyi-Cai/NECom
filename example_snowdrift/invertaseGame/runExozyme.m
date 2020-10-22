% In the case of epsCap = 0, any pairs with vE_1 + vE_2 = vUub / k are Nash
% equilibria (infinitely many).
% But whenever epsCap > 0, no matter how small it is, equal contribution,
% i.e. vE_1 = vE_2 = vUub / 2k becomes the only Nash equilibrium if it is
% able to grow (corresponding to the case epsCap > cost in the linear case in 
% Jeff Gore's original paper about invertase).
% And when the cost is too high or espCap too low, it doesn't grow,
% cooresponding to the case in the original paper where cost > epsCap
%
% Change the parameter pair (epsCap, costExo) to
% 1. (0, 1)     [vE_1 + vE_2 = vUub / k are all Nash equilibria (infinitely many)]
% 2. (0.01, 1)  [Unique NE, all cooperating]
% 3. (0, 7)     [cheating has advantage leading to extinction]
% 4. (0.01, 7)  [same as 3, cheating has advantage leading to extinction] 


N = 2;  % N: number of mutants

% epsCap: capture efficiency 
% (advantage for invertase producer, 0 = none, 1 = get all it breaks down)
epsCap = 0;

% cost of invertase
costExo = 7;

% A pseudo-kinetic parameter
% k_effective for invertase divided by the degradation rate constant
% assuming steady state for the investase in the media
k = 10;
atpm = 0;  % atpm: NGAM value
vUub = 1;  % maximum specific glucose uptake rate

fixMutantVE = [0.05; 0.05];
printAll = 0;
findAllNE = true;
params = struct('OutputFlag', 0);

% find population steady state or find NE at given abundance
findPSS = false;

%X0 = repmat(1/N, N, 1);  % relative abundance
x0 = .5;%rand(N - 1, 1);
X = zeros(numel(x0) + 1, 1);
for j = 1:numel(x0)
    X(j) = prod(x0(1:(j - 1))) * (1 - x0(j));
end
X(end) = prod(x0);
X = flip(X);

% outer-level objective function (weight max for biomass)
outObj = 1 ./ X(:); 
outObj2 = [1 / X(1); zeros(N - 1, 1)];

if ~findPSS
    % Make payoff matrix
    vEvaluesForMatrix = 0:0.01:0.1;%0:0.002:0.1;%[0:0.025:0.075, 0.092:0.002:0.102]';
    nStep = numel(vEvaluesForMatrix);
    
    [gr, vErange, combinations, info0] = exozymeLinearMatrix(X, epsCap, costExo, k, atpm, vUub, vEvaluesForMatrix, params);
    gr1=full(sparse(combinations(:,1),combinations(:,2),gr(:,1),nStep,nStep));
    gr2=full(sparse(combinations(:,1),combinations(:,2),gr(:,2),nStep,nStep));
    [sol, info] = exozymeLinear(X, epsCap, costExo, k, atpm, vUub, 0, outObj, findAllNE, params);
    
%     [sol1, info1] = exozymeLinear(X, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
%     
%     [sol2, info2] = exozymeLinear(X, epsCap, costExo, k, atpm, vUub, 0, outObj, params);
    %% print the pay off matrix
    fprintf('eps = %.2f,  costExo = %.2f,  k = %.2f,  atpm = %.2f,  vUub = %.2f. Linear\n', ...
        epsCap, costExo, k, atpm, vUub);
    fprintf('Payoff Matrix:\n')
    fprintf('Sp1\\Sp2');
    for j = 1:numel(vEvaluesForMatrix)
        fprintf(' %12.04f ', vEvaluesForMatrix(j))
    end
    for j = 1:numel(vEvaluesForMatrix)
        fprintf('\n%6.04f ', vEvaluesForMatrix(j));
        for k2 = 1:numel(vEvaluesForMatrix)
            fprintf(' (%5.03f,%5.03f)', gr1(j, k2), gr2(j, k2));
        end
    end
    fprintf('\n')
    
    LP = info.LP;
    colname = info.colname;
    rowname = info.rowname;
    c = info.c;
    v = info.v;
    fprintf('\nNo mutant defined, max sum of growth rate:\n')
    printSol(X, N, sol, LP, colname, rowname, c, v, printAll)
%     fprintf('\nMutants defined, max sum of growth rate:\n')
%     printSol(epsCap, costExo, k, atpm, vUub, X, N, sol1, info1.LP, colname, rowname, c, v, printAll)
%     fprintf('\nNo mutant defined, max sum of growth rate:\n')
%     printSol(epsCap, costExo, k, atpm, vUub, X, N, sol2, info2.LP, colname, rowname, c, v, printAll)
    
    
else

    %%
    if N > 1
        %X0 = repmat(1/N, N, 1);  % relative abundance
        x0 = rand(N - 1, 1);
        
        if N == 2
            [x, fval, flag, out] = fzero(@(x1) populationSSN2(x1, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params), [0.001, 0.999], optimset('Display','iter'));
            [ss, X, sol, info] = populationSSN2(x, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
            [ss2, X2, sol2, info2] = populationSSN2(0.01, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
            [ss3, X3, sol3, info3] = populationSSN2(0.99, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
        else
            % [ss0, X0, sol0, info0] = populationSSN(x0, epsCap, costExo, k, r, vUub);
            x2 = simulannealbnd(@(z) populationSSN(z, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params), x0, 0.001*ones(N-1, 1), 0.999*ones(N-1, 1));
            [ss, X, sol, info] = populationSSN(x2, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
        end
    else
        [sol, info] = exozymeLinear(1, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
    end
    LP = info.LP;
    colname = info.colname;
    rowname = info.rowname;
    c = info.c;
    v = info.v;
    fprintf('\nNo mutant defined, max sum of growth rate:\n')
    printSol(epsCap, costExo, k, atpm, vUub, X, N, sol, LP, colname, rowname, c, v, printAll)
    printSol(epsCap, costExo, k, atpm, vUub, [0.01; 0.99], N, sol2, info2.LP, colname, rowname, c, v, printAll)
    printSol(epsCap, costExo, k, atpm, vUub, [0.99; 0.01], N, sol3, info3.LP, colname, rowname, c, v, printAll)
end



function printSol(X, N, sol, LP, colname, rowname, c, v, printAll)
if numel(sol) > 1
    for j = 1:numel(sol)
        printSol(X, N, sol(j), LP, colname, rowname, c, v, printAll)
    end
    return
end
if sol.stat == 1
    fprintf('X = (')
    for j = 1:N
        fprintf('%.4f', X(j));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\nmu = (');
    for j = 1:N
        fprintf('%.4f', sol.full(v.vB(j)) ./ X(j));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\nvE/X = (');
    for j = 1:N
        fprintf('%.4f', sol.full(v.vE(j)) ./ X(j));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\nvE/vU = (');
    for j = 1:N
        fprintf('%.4f', sol.full(v.vE(j)) ./ sol.full(v.vU(j)));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\n')
    
    if printAll
        printModel(LP, 'all', 'col', colname, 'row', rowname,'sol',sol)
    else
        printModel(LP, [], 'obj', 'col', colname, 'row', rowname,'sol',sol)
        printModel(LP, [c.mb,c.eflux], 'con', 'col', colname, 'row', rowname,'sol',sol)
    end
    
else
    fprintf('Seems no optimal solution found!\n')
end
end

function [ss, X, sol, info] = populationSSN(x, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params)
X = zeros(numel(x) + 1, 1);
for j = 1:numel(x)
    X(j) = prod(x(1:(j - 1))) * (1 - x(j));
end
X(end) = prod(x);
[ss, sol, info] = populationSS(X, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params);
if isfield(sol, 'objval') && ~isempty(sol.objval)
    if abs(sol.objval) < 1e-5
        ss = 100;
    end
end
end


function [ss, X, sol, info] = populationSSN2(x1, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params)
X = [x1, 1 - x1];
[sol, info] = exozymeLinear(X, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, 0, params);
v = info.v;
if ~strcmp(sol.status, 'OPTIMAL')
    ss = 999;
else
    %     ss = (sol.full(v.vB)' * sol.full(v.vB)) * (X(:)' * X(:)) - (sol.full(v.vB)' * X(:)) ^ 2;
    ss = sol.full(v.vB) ./ X(:);
    if any(abs(ss) < 1e-4)
        ss = 999;
    else
        ss = ss(2) - ss(1);
    end
end
end

function [ss, sol, info] = populationSS(X, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, params)
[sol, info] = exozymeLinear(X, epsCap, costExo, k, atpm, vUub, fixMutantVE, outObj, 0, params);
v = info.v;
if ~strcmp(sol.status, 'OPTIMAL')
    ss = 999;
else
    %     ss = (sol.full(v.vB)' * sol.full(v.vB)) * (X(:)' * X(:)) - (sol.full(v.vB)' * X(:)) ^ 2;
    ss = norm((sol.full(v.vB) ./ X(:)) - mean(sol.full(v.vB) ./ X(:)));
end
end
