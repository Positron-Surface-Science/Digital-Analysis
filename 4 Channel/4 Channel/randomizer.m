for i=1:numel(p4s(1,:))
    x = rand;
    
    if x < 0.5
        temp = p4s(:,i);
        p4s(:,i) = p4s(:,numel(p4s(1,:))-i);
         p4s(:,numel(p4s(1,:))-i) = temp;
        
    end
    
end