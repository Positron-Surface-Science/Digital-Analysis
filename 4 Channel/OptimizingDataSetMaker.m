x = linspace(1,438,438)';
nt = [];
targets = [];
%netMatrix = net(tN(:,1:27000));
no = 169

for j=no:no
    
   % if exact(j) == 1
        
        nt = horzcat(nt,tN(:,netMatrix(j,:) == 1));
        
        %tstest = zeros(2,numel(tN(1,netMatrix(j,:) == 1))) + ELETMatrix(j,1:2)';
        
        out = net(ntTest);
        
        %best = mean(out,2);
        
        %ELETMatrixBase(j,:) = best';
        
        %targets = horzcat(targets,ts);
        
        
        plot(ax1,x,nt)
        xticks(ax1,0:10:400);
        grid(ax1,'on')
       %{
        beep
        j
        yesno = input('Good Traces? (y/n): ','s');
        
        if yesno == 'n'
            exact(j) = 0;
            
        end
        %}
        %{
        if yesno == 'y'
        
            template = mean(nt,2);
            template = template/max(template);
        derivative = diff(template);
        ds = smooth(derivative,13,'sgolay');
        d2 = smooth(diff(ds),13,'sgolay');
        plot(ax2,x(1:400),derivative)
        xticks(ax2,0:10:400);
        grid(ax2,'on')
        
        range1 = round(range1 + 0.25);%input('Range 1: ');
        range2 = round(range2 + 0.25);%input('Range 2: ');
        
        
            %residuals = (nt - template).^2;
            
            %[minimum,minIndex] = max(mean(residuals(range1:range2,:),2));
            
            [maxDiff,maxdIndex] = min(d2(range1:range2));
            [minDiff,mdIndex] = min(ds(range1:range2));
            
            mdIndex = mdIndex + range1;
            maxdIndex = maxdIndex + range1;
            minIndex = minIndex + range1;
            
            %minMean = round(minIndex);
            %minAve = mean(residuals(range1:range2,:));
            
            percentage(j,:) = [template(maxdIndex)*100 template(mdIndex)*100];
            'Percentage: '
            percentage(j,:)
            
        end
        %}
    %end
end