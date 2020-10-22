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
% 1. (0, 10)     [vE_1 + vE_2 = vUub / k are all Nash equilibria (infinitely many), all cooperating has higher sum of growth]
% 2. (0.01, 10)  [Unique NE, all cooperating]
% 3. (0.01, 30)  [similar to case 1, multiple NE]
% 4. (0, 50)     [cheating has advantage leading to extinction]
% 5. (0.01, 50)  [same as case 4, cheating has advantage leading to extinction] 

N = 2;  % N: number of mutants

% epsCap: capture efficiency 
% (advantage for invertase producer, 0 = none, 1 = get all it breaks down)
epsCap = 0.;

% cost of invertase
costInv = 50;

% quadratic cost of invertase
qc = 2;

% A pseudo-kinetic parameter
% k_effective for invertase divided by the degradation rate constant
% assuming steady state for the investase in the media
k = 100;
r = 10;  % atpm: NGAM value
vUub = 10;  % maximum specific glucose uptake rate

fixMutantVE = false; %[0; 0.05];
findAllNE = true;
printAll = 0;
params = struct('OutputFlag', 0);

% find population steady state or find NE at given abundance
findPSS = false;

% Create an initial relative abundance profile for N members using N - 1 numbers
x0 = 0.5; % rand(N - 1, 1);
assert(numel(x0) == N - 1)
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
    vEvaluesForMatrix = 0:(vUub / k / 10):(vUub / k);%0:0.002:0.1;%[0:0.025:0.075, 0.092:0.002:0.102]';
    nStep = numel(vEvaluesForMatrix);
    
    [gr, vErange, combinations, info0] = exozymeQuadMatrix(X, epsCap, costInv, k, r, vUub, qc, vEvaluesForMatrix, params);
    gr1=full(sparse(combinations(:,1),combinations(:,2),gr(:,1),nStep,nStep));
    gr2=full(sparse(combinations(:,1),combinations(:,2),gr(:,2),nStep,nStep));
    [sol, info] = exozymeQuad(X, epsCap, costInv, k, r, vUub, qc, 0, outObj, findAllNE, params);
    
%     [sol1, info1] = exozymeQuad(X, epsCap, costInv, k, r, vUub, qc, [0; NaN],outObj2, params);
%     
%     [sol2, info2] = exozymeQuad(X, epsCap, costInv, k, r, vUub, qc, 0,outObj2, params);
    %% print the pay off matrix
    fprintf('eps = %.2f,  costInv = %.2f,  k = %.2f,  r = %.2f,  vUub = %.2f. Quadratic\n', ...
        epsCap, costInv, k, r, vUub);
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
%     fprintf('\nMutant Sp1 defined, max sp1 growth rate:\n')
%     printSol(epsCap, costInv, k, r, vUub, X, N, sol1, info1.LP, colname, rowname, c, v, printAll)
%     fprintf('\nNo mutant defined, max sp1 growth rate:\n')
%     printSol(epsCap, costInv, k, r, vUub, X, N, sol2, info2.LP, colname, rowname, c, v, printAll)
    
    
else

    %%
    if N > 1
        %X0 = repmat(1/N, N, 1);  % relative abundance
        x0 = rand(N - 1, 1);
        
        if N == 2
            [x, fval, flag, out] = fzero(@(x1) populationSSN2(x1, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params), [0.001, 0.999], optimset('Display','iter'));
            [ss, X, sol, info] = populationSSN2(x, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params);
            [ss2, X2, sol2, info2] = populationSSN2(0.01, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params);
            [ss3, X3, sol3, info3] = populationSSN2(0.99, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params);
        else
            % [ss0, X0, sol0, info0] = populationSSN(x0, epsCap, costInv, k, r, vUub);
            x2 = simulannealbnd(@(z) populationSSN(z, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params), x0, 0.001*ones(N-1, 1), 0.999*ones(N-1, 1));
            [ss, X, sol, info] = populationSSN(x2, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params);
        end
    else
        [sol, info] = exozymeQuad(1, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, 0, params);
    end
    LP = info.LP;
    colname = info.colname;
    rowname = info.rowname;
    c = info.c;
    v = info.v;
    fprintf('\nNo mutant defined, max sum of growth rate:\n')
    printSol(X, N, sol, LP, colname, rowname, c, v, printAll)
    printSol([0.01; 0.99], N, sol2, info2.LP, colname, rowname, c, v, printAll)
    printSol([0.99; 0.01], N, sol3, info3.LP, colname, rowname, c, v, printAll)
end



function printSol(X, N, sol, LP, colname, rowname, c, v, printAll)
if numel(sol) > 1
    for j = 1:numel(sol)
        printSol(X, N, sol(j), LP, colname, rowname, c, v, printAll)
    end
    return
end
if strcmp(sol.status, 'OPTIMAL')
    fprintf('X = (')
    for j = 1:N
        fprintf('%.4f', X(j));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\nmu = (');
    for j = 1:N
        fprintf('%.4f', sol.x(v.vB(j)) ./ X(j));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\nvE/X = (');
    for j = 1:N
        fprintf('%.4f', sol.x(v.vE(j)) ./ X(j));
        if j ~= N
            fprintf(',');
        end
    end
    fprintf(')\nvE/vU = (');
    for j = 1:N
        fprintf('%.4f', sol.x(v.vE(j)) ./ sol.x(v.vU(j)));
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

function [ss, X, sol, info] = populationSSN(x, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params)
X = zeros(numel(x) + 1, 1);
for j = 1:numel(x)
    X(j) = prod(x(1:(j - 1))) * (1 - x(j));
end
X(end) = prod(x);
[ss, sol, info] = populationSS(X, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params);
if isfield(sol, 'objval') && ~isempty(sol.objval)
    if abs(sol.objval) < 1e-5
        ss = 100;
    end
end
end


function [ss, X, sol, info] = populationSSN2(x1, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params)
X = [x1, 1 - x1];
[sol, info] = exozymeQuad(X, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, 0, params);
v = info.v;
if ~strcmp(sol.status, 'OPTIMAL')
    ss = 999;
else
    %     ss = (sol.x(v.vB)' * sol.x(v.vB)) * (X(:)' * X(:)) - (sol.x(v.vB)' * X(:)) ^ 2;
    ss = sol.x(v.vB) ./ X(:);
    if any(abs(ss) < 1e-4)
        ss = 999;
    else
        ss = ss(2) - ss(1);
    end
end
end

function [ss, sol, info] = populationSS(X, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, params)
[sol, info] = exozymeQuad(X, epsCap, costInv, k, r, vUub, qc, fixMutantVE, outObj, 0, params);
v = info.v;
if ~strcmp(sol.status, 'OPTIMAL')
    ss = 999;
else
    %     ss = (sol.x(v.vB)' * sol.x(v.vB)) * (X(:)' * X(:)) - (sol.x(v.vB)' * X(:)) ^ 2;
    ss = norm((sol.x(v.vB) ./ X(:)) - mean(sol.x(v.vB) ./ X(:)));
end
end
