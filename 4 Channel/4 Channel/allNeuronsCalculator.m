xToF = linspace(-0.25E-6,0.25E-6,1024)';

for n=1:225
    
    if exact(n) == 1
        
        for i=1:220
            
            if allNeurons.eletParameters(1,i) ~= 0
                
                allNeurons.neurons{n}(allNeurons.neurons{n}(:,i) == 0,i) = NaN;
                
                hToF = histcounts(allNeurons.neurons{n}(:,i),'NumBins',1024,'BinLimits',[-0.25E-6 0.25E-6])';
                
                observe(i,n) = numel(allNeurons.neurons{n}(isnan(allNeurons.neurons{n}(:,i)) == 0,i));
                
                if max(hToF) < 200
                    
                    try
                        gaussFit = fit(xToF,hToF,'gauss1');
                        
                        cValues = coeffvalues(gaussFit);
                        
                        fwhm(i,n) = (2*sqrt(log(2))*cValues(3))/sqrt(2);
                        amplitude(i,n) = cValues(1);
                        x0(i,n) = cValues(2);
                        
                    catch
                        fwhm(i,n) = 1;
                        amplitude(i,n) = 0;
                        
                    end
                    
                end
                
            end
            
        end
        
        [a,b] = min(fwhm(observe(:,n) >= 0.975*max(observe(:,n)),n));
        
        index = find(fwhm(:,n) <= 1.025*a);
        
        [~,minFwhmIndex] = min(fwhm(index,n));
        
        smallestIndex = index(minFwhmIndex);
        
        if a ~= 0
            best(:,n) = vertcat(allNeurons.eletParameters(:,smallestIndex),x0(smallestIndex,n));
            
        else
            best(:,n) = [0 0 0];
            
        end
        
    else
        best(:,n) = [0 0 0];
        
    end
    
    n
    
end