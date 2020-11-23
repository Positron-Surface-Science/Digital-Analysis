range1 = 250;
range2 = 1750;

specRange1 = 1;
specRange2 = 75;

t1 = 15;
t2 = 5;
t3 = 250;

nPulsesPerPeak = 4250;

fs = 2.5E9;

data = [];
response = [];

for i=1:nPulsesPerPeak
    
    spec = spectrogram(p1(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
        data{i}(:,:) = rispec*1E3;
        
    
        response{i} = single(p1(range1:range2,i)'*1E3);
        
    
    
    
    spec = spectrogram(p2(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
        data{i + nPulsesPerPeak}(:,:) = rispec*1E3;
        
    
        response{i + nPulsesPerPeak} = single(p2(range1:range2,i)'*1E3);
        
    
    
    spec = spectrogram(p3(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
        data{i + 2*nPulsesPerPeak}(:,:) = rispec*1E3;
        
    
        response{i + 2*nPulsesPerPeak} = single(p3(range1:range2,i)'*1E3);
        
   
    
    spec = spectrogram(p4(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
        data{i + 3*nPulsesPerPeak}(:,:) = rispec*1E3;
        
    
        response{i + 3*nPulsesPerPeak} = single(p4(range1:range2,i)'*1E3);
        
    
    
    spec = spectrogram(p5(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
        data{i + 4*nPulsesPerPeak}(:,:) = rispec*1E3;
        
    
        response{i + 4*nPulsesPerPeak} = single(p5(range1:range2,i)'*1E3);
        
    
    
end

for i=1:(nPulsesPerPeak)*5
    response{i} = response{i}';
    
end
