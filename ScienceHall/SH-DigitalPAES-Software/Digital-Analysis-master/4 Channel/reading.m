fid = fopen('E:\DRS4\Data\Ge\May\Ba-133_Eu152\waveform.dat');

for i=1:10000
[amp(i),shapedPulse(i,:)] = DRS4Read(fid,i);
i
end