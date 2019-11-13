
f1 = figure;
ax1 = axes(f1);
f2 = figure;
ax2 = axes(f2);

%x = linspace(1,401,401);

%percentage = zeros(144,2);
range1 = 160;
range2 = 180;

for j=1:144
    
    nt=[];
    
    if exact(j) == 1 && percentage2(j) ~= 0 %&& ...
            %(percentage(j) < 0.01 || percentage(j) >= 0.95)
            
        netMatrix = net(treN);
        
        nt = treN(:,netMatrix(j,:) == 1);
        
        plot(ax1,x,nt)
        xticks(ax1,0:10:400);
        grid(ax1,'on')
        beep
        j
        %yesno = input('Good Traces? (y/n): ','s');
        
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
    end
end