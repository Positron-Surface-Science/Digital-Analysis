range1 = 200;
range2 = 2500;

specRange1 = 10;
specRange2 = 90;

fs = 2.5E9;

for i=1:4000
    
    spec = spectrogram(p1(range1:range2,1), 5, 2, 250, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= 2000
        data{i}(:,:) = rispec*1E4;
        
    else
        response{i - 2000} = p1(range1:range2,i)'*1E4;
        
    end
    
end

for i=1:4000
    
    spec = spectrogram(p2(range1:range2,1), 5, 2, 250, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= 2000
        data{i + 2000}(:,:) = rispec*1E4;
        
    else
        response{i} = p2(range1:range2,i)'*1E4;
        
    end
    
end


for i=1:4000
    
    spec = spectrogram(p3(range1:range2,1), 5, 2, 250, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= 2000
        data{i + 4000}(:,:) = rispec*1E4;
        
    else
        response{i + 2000} = p3(range1:range2,i)'*1E4;
        
    end
    
end

for i=1:4000
    
    spec = spectrogram(p4(range1:range2,1), 5, 2, 250, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= 2000
        data{i + 6000}(:,:) = rispec*1E4;
        
    else
        response{i + 4000} = p4(range1:range2,i)'*1E4;
        
    end
    
end

for i=1:4000
    
    spec = spectrogram(p5(range1:range2,1), 5, 2, 250, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= 2000
        data{i + 8000}(:,:) = rispec*1E4;
        
    else
        response{i + 6000} = p5(range1:range2,i)'*1E4;
        
    end
    
end

for i=1:10000
    response{i} = response{i}';
    
end
