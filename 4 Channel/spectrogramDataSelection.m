

for i=1:2000
    
    spec = spectrogram(p1(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 1000
        data{i}(:,:) = rispec;
        
    else
        response{i - 1000}(:,:) = rispec;
        
    end
    
end

for i=1:2000
    
    spec = spectrogram(p2(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 1000
        data{i + 1000}(:,:) = rispec;
        
    else
        response{i}(:,:) = rispec;
        
    end
    
end


for i=1:2000
    
    spec = spectrogram(p3(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 1000
        data{i + 2000}(:,:) = rispec;
        
    else
        response{i + 1000}(:,:) = rispec;
        
    end
    
end

for i=1:2000
    
    spec = spectrogram(p4(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 1000
        data{i + 3000}(:,:) = rispec;
        
    else
        response{i + 2000}(:,:) = rispec;
        
    end
    
end

for i=1:2000
    
    spec = spectrogram(p5(350:1750,i), 100, 'yaxis');
    
    rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
    
    if i <= 1000
        data{i + 4000}(:,:) = rispec;
        
    else
        response{i + 3000}(:,:) = rispec;
        
    end
    
end
