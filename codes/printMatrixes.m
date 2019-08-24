% Tradeoff matrix print

function printMatrixes(theresults,theVarRxns,EX_Rflux,strsize,filename)
% Function: write the payoff matrices for the state transition map
% Jingyi Cai 2019-5

if (nargin < 4)
    strsize = [4,4];
end
if (nargin < 5)
    filename = 'TradeOffMatrixes.xls';
end

matrixSize(1)=length(EX_Rflux(1,:));
matrixSize(2)=length(EX_Rflux(2,:));
% get the tradeoff matrix

TM=theresults.Tradoffs;

h=waitbar(0,'writing...');
steps=matrixSize(1)*matrixSize(2) ;
stepcount=0;

for i=1:matrixSize(1)
    for j=1:matrixSize(2) 
        stepcount=stepcount+1;
        matrixSt{i,j}=[(i-1)*(strsize(1)+2)+3,(j-1)*(strsize(1)+1)+3];
        theRange=conv2xlscell(matrixSt{i,j}(1), matrixSt{i,j}(2));
        if ~isempty(TM{i,j})
            theMatrx=cellfun(@(x) mat2str(x),TM{i,j}, 'UniformOutput', false);
        else
            theMatrx=[];
        end
        xlswrite(filename,theMatrx,'TradeoffM',theRange);
        if i==1
            theRange2=conv2xlscell(matrixSt{i,j}(1)-2, matrixSt{i,j}(2)+2);
            xlswrite(filename,EX_Rflux(2,j),'TradeoffM',theRange2);
            theRange2_=conv2xlscell(matrixSt{i,j}(1)-1, matrixSt{i,j}(2));
            xlswrite(filename,{'1,1','1,0','0,1','0,0'},'TradeoffM',theRange2_);
        end
        if j==1
            theRange3=conv2xlscell(matrixSt{i,j}(1)+2, matrixSt{i,j}(2)-2);
            xlswrite(filename,EX_Rflux(1,i),'TradeoffM',theRange3);
            theRange3_=conv2xlscell(matrixSt{i,j}(1), matrixSt{i,j}(2)-1);
            xlswrite(filename,{'1,1';'1,0';'0,1';'0,0'},'TradeoffM',theRange3_);
        end   
        waitbar(stepcount/steps);
    end
end

xlswrite(filename,theVarRxns(1),'TradeoffM','A2');
xlswrite(filename,theVarRxns(2),'TradeoffM','B1');
close(h)
