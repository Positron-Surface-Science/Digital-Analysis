range1 = 250;
range2 = 1750;

specRange1 = 1;
specRange2 = 150;

t1 = 10;
t2 = 2;
t3 = 350;

nPulsesPerPeak = 4250;

fs = 2.5E9;

data = [];
response = [];

for i=1:nPulsesPerPeak
    
    spec = spectrogram(p1(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    %if i <= nPulsesPerPeak/2
        data{i}(:,:) = rispec*1E3;
        
    %else
        response{i - nPulsesPerPeak/2} = single(p1(range1:range2,i)'*1E3);
        
    %end
    
    
    spec = spectrogram(p2(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    %if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak/2}(:,:) = rispec*1E3;
        
    %else
        response{i} = single(p2(range1:range2,i)'*1E3);
        
    %end
    
    spec = spectrogram(p3(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    %if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak}(:,:) = rispec*1E3;
        
    %else
        response{i + nPulsesPerPeak/2} = single(p3(range1:range2,i)'*1E3);
        
    %end
    
    spec = spectrogram(p4(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    %if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak/2 + nPulsesPerPeak}(:,:) = rispec*1E3;
        
    %else
        response{i + nPulsesPerPeak} = single(p4(range1:range2,i)'*1E3);
        
    %end
    
    spec = spectrogram(p5(range1:range2,1), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    %if i <= nPulsesPerPeak/2
        data{i + nPulsesPerPeak*2}(:,:) = rispec*1E3;
        
    %else
        response{i + nPulsesPerPeak/2 + nPulsesPerPeak} = single(p5(range1:range2,i)'*1E3);
        
    %end
    
end

for i=1:(nPulsesPerPeak/2)*5
    response{i} = response{i}';
    
end
