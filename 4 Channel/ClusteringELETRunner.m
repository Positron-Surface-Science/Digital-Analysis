
%f1 = figure;
%ax1 = axes(f1);
nt.x = linspace(1,438+18626,438+18626);
nt.desc.Ts = 1/2.5E9;
%netMatrix = net(tN(:,1:75000));
DV.exact = exact;
%best = zeros(225,2);
maxima = [];

for j=1:225
    
    nt.y = [];
    ELETMatrix(j,:) = [1 3];
    
    if exact(j) == 1
        
        all = tN(:,netMatrix(j,:) == 1);
        numberRun = 1;
        
        for u=1:95
            
            for s=1:100
                
                for n=1:numel(all(1,:))
                    
                    nt.y = vertcat(zeros(18626,1),all(:,n))';
                    
                    [crossTimeOut(n),noGood,numPeaks,~] = paes(0,3,1,0,...
                        2,1,nt,0,1E9,ELETMatrix,DV,j,j,1,1);
                    
                    %{
                    [crossTimeOut,noGood(c,i),numPeaks,~] = paes(timingTypeIn(c),TypeIn(c),c,multiStopIn,...
                    multiStopConditionsIn,numPeaks,pulse{c}(s),RTFLIn,RTFUpIn,ELETMatrix, ...
                    DiscriminationVector,neuron,activeNeuron,ANN,ggc);
                    %}
                end
                
                hToF = histcounts(crossTimeOut,'NumBins',1024,'BinLimits', ...
                        [1 20000])';
                    
                    xToF = linspace(1,20000,1024)';
                    
                    max2 = max(hToF);
                    
                    if numberRun > 1 && max(hToF) > max(maxima)
                        
                        best(j,:) = ELETMatrix(j,:)';
                        
                    end
                    
                    maxima(numberRun) = max2;
                    
                    ELETMatrix(j,2) = ELETMatrix(j,2) + 1;
                    ELETMatrix(j,:);
                    
                    
                    
                    numberRun = numberRun + 1;
                    
            end
            
            ELETMatrix(j,1) = ELETMatrix(j,1) + 1;
            ELETMatrix(j,2) = ELETMatrix(j,1) + 1;
            
        end
        
    end
    
    j
    
end