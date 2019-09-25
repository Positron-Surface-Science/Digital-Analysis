for c=1:100
    
    [m,n] = size(geHistogram{c});
    
    input = horzcat(input,geHistogram{c}(:,1:floor(n/2)));
    output = horzcat(output,geHistogram{c}(:,floor(n/2):floor(n/2)*2-1));
    
    
    c
    
end


