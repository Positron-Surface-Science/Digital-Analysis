%for j=1:100
    
    nt=[];
    
    %if output2(j) >= 2500
        
        for i=1:1000
            
            compares = ClusteringNetwork.net(treN(:,i));
            
            if compares(14) == 1
                
                nt = horzcat(nt,treN(:,i));
                
            end
            
        end
        
        figure
        plot(nt)
        
    %end
%end