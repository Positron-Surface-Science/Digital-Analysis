b1 = mean(p1s(200:400,:));
p1s = p1s - b1;
b2 = mean(p2s(200:400,:));
p2s = p2s - b2;
b3 = mean(p3s(200:400,:));
p3s = p3s - b3;
b4 = mean(p4s(200:400,:));
p4s = p4s - b4;
b5 = mean(p5s(200:400,:));
p5s = p5s - b5;

range1 = 501;
range2 = 950;

specRange1 = 1;
specRange2 = 176;

t1 = 5;
t2 = 2;
t3 = 350;

nPulsesPerPeak = 16000;

fs = 2.5E9;

data = [];
response = [];

pAverage = ps1Mean; %mean(mean(p1s(range1:range2,i)));

pStd = ps1Std; std(std(p1s(range1:range2,i)));

for i=1:nPulsesPerPeak
    
    spec = spectrogram(p1s(range1:range2,i), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    average = mean(mean(rispec));
    
    std = std(std(rispec));
    
    data{i}(:,:) = (rispec - average)/std;
    
    
    response{i} = single((p1s(range1:range2,i)' - pAverage)/pStd);
    
    
    
    
    spec = spectrogram(p2s(range1:range2,i), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
    average = mean(mean(rispec));
    
    std = std(std(rispec));(rispec - average)/std;
    
    data{i + nPulsesPerPeak}(:,:) = (rispec - average)/std;
    
    response{i + nPulsesPerPeak} = single((p2s(range1:range2,i)' - pAverage)/pStd);
    
    
    
    spec = spectrogram(p3s(range1:range2,i), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    
    average = mean(mean(rispec));
    
    std = std(std(rispec));
    
    data{i + 2*nPulsesPerPeak}(:,:) = (rispec - average)/std;
    
    response{i + 2*nPulsesPerPeak} = single((p3s(range1:range2,i)' - pAverage)/pStd);
    
    
    
    spec = spectrogram(p4s(range1:range2,i), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
     average = mean(mean(rispec));
    
    std = std(std(rispec));
    
    data{i + 3*nPulsesPerPeak}(:,:) = (rispec - average)/std;
    
    
    response{i + 3*nPulsesPerPeak} = single((p4s(range1:range2,i)' - pAverage)/pStd);
    
    
    
    
    spec = spectrogram(p5s(range1:range2,i), t1, t2, t3, fs, 'yaxis');
    
    rispec = horzcat(single(real(spec(specRange1:specRange2,:))), single(imag(spec(specRange1:specRange2,:))));
    
    average = mean(mean(rispec));
    
    std = std(std(rispec));
    
    data{i + 4*nPulsesPerPeak}(:,:) = (rispec - average)/std;
    
    response{i + 4*nPulsesPerPeak} = single((p5s(range1:range2,i)' - pAverage)/pStd);
    
    
    
end

for i=1:(nPulsesPerPeak)*5
    response{i} = response{i}';
    
end
