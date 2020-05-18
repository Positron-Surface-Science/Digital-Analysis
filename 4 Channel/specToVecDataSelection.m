

for i=1:4000
    
    spec = spectrogram(p1(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 2000
        data{i}(:,:) = rispec;
        
    else
        response{i - 2000} = p1(350:1750,i)';
        
    end
    
end

for i=1:4000
    
    spec = spectrogram(p2(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 2000
        data{i + 2000}(:,:) = rispec;
        
    else
        response{i} = p2(350:1750,i)';
        
    end
    
end


for i=1:4000
    
    spec = spectrogram(p3(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 2000
        data{i + 4000}(:,:) = rispec;
        
    else
        response{i + 2000} = p3(350:1750,i)';
        
    end
    
end

for i=1:4000
    
    spec = spectrogram(p4(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 2000
        data{i + 6000}(:,:) = rispec;
        
    else
        response{i + 4000} = p4(350:1750,i)';
        
    end
    
end

for i=1:4000
    
    spec = spectrogram(p5(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 2000
        data{i + 8000}(:,:) = rispec;
        
    else
        response{i + 6000} = p5(350:1750,i)';
        
    end
    
end

for i=1:10000
    response{i} = response{i}';
    
end
