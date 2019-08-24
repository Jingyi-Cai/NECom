function theresults=StateTransitionMap(themodel,parameters2,ctrlinds,EX_Rflux)


%{
Function: prepare all payoff matrices for the state transition Map
Input: 
    themodel: microbial community model struct
    parameters2: parameters
    ctrlinds: indices of the uptake reactions to be controlled to capture state for StateTransitionMap. 
    EX_Rflux: controlling fluxes for the controlled uptake reactions.
Output: 
    theresults:contains all payoff matrice(FSHaP) for the state transition map 
Jingyi Cai 2019-08
%}



EX_Rflux1=EX_Rflux{1,1};
EX_Rflux2=EX_Rflux{1,2};
Tradoffs=[];
SpRxnExFluxes=[];

h=waitbar(0,'computing...');
steps=length(EX_Rflux1)*length(EX_Rflux2);
stepcount=0;
for i=1:length(EX_Rflux1)
    for j=1:length(EX_Rflux2)
        stepcount=stepcount+1;
        themodel_ed=themodel;
        themodel_ed.lb(ctrlinds)=[EX_Rflux1(i),EX_Rflux2(j)];
        themodel_ed.ub(ctrlinds)=[EX_Rflux1(i),EX_Rflux2(j)];    
        % capture a state using joint-FBA
        solfba1=fba(themodel_ed); % Flux Balance Analysis
        x1=solfba1.x;
        % then construct the payoff matrix(FShaP) at this state. 
        if ~isempty(x1)
            SpRxnExFluxes{i,j}=x1(parameters2.sub_indExSpi);     
            [Tradoffs{i,j},theflagss{i,j}]=...
                FShaP(themodel,parameters2,SpRxnExFluxes{i,j});
        else
            SpRxnExFluxes{i,j}=[];
        end
        waitbar(stepcount/steps);
    end
end
theresults.Tradoffs=Tradoffs;
theresults.theflagss=theflagss;
close(h)




end







