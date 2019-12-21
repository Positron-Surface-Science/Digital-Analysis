
f1 = figure;
ax1 = axes(f1);
f2 = figure;
ax2 = axes(f2);
nt=[];
nt.x = linspace(1,438+18626,438+18626);
nt.desc.Ts = 1/2.5E9;
%netMatrix = net(tN(:,1:75000));
DV.exact = exact;
%best = zeros(225,2);
maxima = [];
base = [];
ELETMatrix = cell(1,1500);
elet = zeros(225,2);
hToFs = [];

for j=85:85
    
    if exact(j) == 1
        
        numberRun = 1;
        allSigma = [];
        maxima = [];
        crossTimeOut = [];
        nt.y = [];
        observations = [];
        allC1 = [];
        allC3 = [];
        
        all = tN(:,netMatrix(j,:) == 1);
        
        base(j) = mean(mean(all(1:10,:)));
        
        low = (round(base(j)*100) + 1);
        high = (round(base(j)*100) + 3);
        
        ELETMatrix{numberRun} = [low high];
        best(j,1:2) = ELETMatrix{numberRun}';
        
        for u=1:40
            
            if ELETMatrix{numberRun}(1) > 100
                break;
                
            end
            
            for s=u+1:40
                
                if ELETMatrix{numberRun}(2) > 100
                    break;
                    
                end
                
                try
                    full = numel(all(1,1:100));
                    
                catch
                    full = numel(all(1,:));
                    
                end
                
                for n=1:full
                    
                    nt.y = vertcat(zeros(18626,1),all(:,n))';
                    
                    elet(j,:) = ELETMatrix{numberRun};
                    
                    [crossTimeOut(n),noGood,numPeaks,~] = paes(0,3,1,0,...
                        2,1,nt,0,1E9,elet,DV,j,j,1,1);
                    
                    %{
                    [crossTimeOut,noGood(c,i),numPeaks,~] = paes(timingTypeIn(c),TypeIn(c),c,multiStopIn,...
                    multiStopConditionsIn,numPeaks,pulse{c}(s),RTFLIn,RTFUpIn,ELETMatrix, ...
                    DiscriminationVector,neuron,activeNeuron,ANN,ggc);
                        %}
                end
                
                hToF = histcounts(crossTimeOut,'NumBins',3000,'BinLimits', ...
                    [18500 20000])';
                
                try
                    
                %hToF = hToF./max(hToF);
                
                xToF = linspace(18500,20000,3000)';
                
                observe = numel(crossTimeOut(isnan(crossTimeOut) == 0));
                
                gaussFit = fit(xToF,hToF,'gauss1');
                
                cValues = coeffvalues(gaussFit);
                
                max2 = max(hToF);
                
                maxima(numberRun) = max2;
                allC3(numberRun) = cValues(3);
                allC1(numberRun) = cValues(1);
                observations(numberRun) = observe;
                
                hToFs(numberRun,:) = hToF;
                
                numberRun = numberRun + 1;
                
                ELETMatrix{numberRun}(1) =  ELETMatrix{numberRun-1}(1);
                ELETMatrix{numberRun}(2) =  ELETMatrix{numberRun-1}(2) + 1;
                
                catch
                    break;
                    
                end
                
            end
            
            ELETMatrix{numberRun}(1) =  ELETMatrix{numberRun-1}(1) + 1;
            ELETMatrix{numberRun}(2) =  ELETMatrix{numberRun}(1) + 2;
            j
            ELETMatrix{numberRun}
            
        end
        
        [a,b] = min(allC3(allC3 >= 1 & observations == max(observations)));
        location = find(allC3 <= 1.01*a & allC3 >= 1 & observations == max(observations));
        best(j,1:2) = ELETMatrix{location};
        pause(1)
        
    end
    
end