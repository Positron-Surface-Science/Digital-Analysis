xCh = linspace(0,7.5,4096)';


    %observe = [];
    %fwhm = [];
    %index = [];
    %best = [];
    

for n=1:400
    
    %if exact(n) == 1
        
        for i=1:108
            
            if allNeurons.mwdParameters(2,i) ~= 0
                
                allNeurons.mwdNeurons{n}(allNeurons.mwdNeurons{n}(:,i) == 0,i) = NaN;
                
                hCh = histcounts(allNeurons.mwdNeurons{n}(:,i),'NumBins',4096,'BinLimits',[0 7.5])';
                hCh = hCh/max(hCh);
                
                observe(i,n) = numel(allNeurons.mwdNeurons{n}(isnan(allNeurons.mwdNeurons{n}(:,i)) == 0,i));
                
                if max(hCh) < sum(hCh) && observe(i,n) >= 0.75*max(observe(:,n))
                    
                    try
                    
                    [~,second] = max(hCh);
                    
                    if second < observe(i,n)
                        gaussFit = fit(xCh(second-250:second+250),hCh(second-250:second+250),'gauss1');
                    else
                        gaussFit = fit(xCh(1:second+250),hCh(1:second+250),'gauss1');
                    end
                    
                        cValues = coeffvalues(gaussFit);
                        
                        fwhm(i,n) = (2*sqrt(log(2))*cValues(3))/sqrt(2);
                        amplitude(i,n) = cValues(1);
                        x0(i,n) = cValues(2);
                        
                        %{
                        nVector = allNeurons.neurons{n}(isnan(allNeurons.neurons{n}(:,i)) == 0 | ~any(allNeurons.neurons{n}(:,i)),i);
                        r = fitdist(nVector,'Normal');
                        fwhm(i,n) = r.sigma;
                        x0(i,n) = r.mu;
                        %}
                        %{
                        nVector = allNeurons.neurons{n}(isnan(allNeurons.neurons{n}(:,i)) == 0 | ~any(allNeurons.neurons{n}(:,i)),i);
                        [mu,sigma] = normfit(nVector);
                        fwhm(i,n) = sigma;
                        x0(i,n) = mu;
                        %}
                        %fwhm(i,n) = std(nVector);
                        %x0(i,n) = mean(nVector);
                        
                    catch
                        fwhm(i,n) = 1;
                        amplitude(i,n) = 0;
                        x0(i,n) = 0;
                        'error'
                        
                    end
                    
                end
                
            end
            
        end
        
        try
        [a,b] = min(fwhm(observe(:,n) >= 0.75*max(observe(:,n)),n));
        
        index = find(fwhm(:,n) <= 1.05*a & fwhm(:,n) > 0);
        
        [~,minFwhmIndex] = min(fwhm(index,n));
        
        smallestIndex = index(minFwhmIndex);
        
        if a ~= 0
            best(:,n) = vertcat(allNeurons.mwdParameters(:,smallestIndex),x0(smallestIndex,n),fwhm(smallestIndex,n));
            'sdgsdaf'
            
        else
            best(:,n) = [0 0 0 0];
            
        end
        
    %else
     %   best(:,n) = [0 0 0 0];
        
    %end
    
    n
    
        catch
            best(:,n) = [0 0 0 0];
        end
    
end