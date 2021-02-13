range1 = 451;
range2 = 950;

specRange1 = 1;
specRange2 = 75;

t1 = 15;
t2 = 10;
t3 = 200;

nPulsesPerPeak = 15000;

fs = 2.5E9;

data = [];
response = [];

for i=1:nPulsesPerPeak
    
    spec = spectrogram(p1s(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= nPulsesPerPeak/2
        data{i}(:,:) = rispec*1E3;
        
    else
        response{i - nPulsesPerPeak/2} = single(p1s(range1:range2,i)'*1E3);
        
    end
    
    
    spec = spectrogram(p2s(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak/2}(:,:) = rispec*1E3;
        
    else
        response{i} = single(p2s(range1:range2,i)'*1E3);
        
    end
    
    spec = spectrogram(p3s(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak}(:,:) = rispec*1E3;
        
    else
        response{i + nPulsesPerPeak/2} = single(p3s(range1:range2,i)'*1E3);
        
    end
    
    spec = spectrogram(p4s(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak/2 + nPulsesPerPeak}(:,:) = rispec*1E3;
        
    else
        response{i + nPulsesPerPeak} = single(p4s(range1:range2,i)'*1E3);
        
    end
    
    spec = spectrogram(p5s(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak*2}(:,:) = rispec*1E3;
        
    else
        response{i + nPulsesPerPeak/2 + nPulsesPerPeak} = single(p5s(range1:range2,i)'*1E3);
        
    end
    
end

for i=1:(nPulsesPerPeak/2)*5
    response{i} = response{i}';
    
end
