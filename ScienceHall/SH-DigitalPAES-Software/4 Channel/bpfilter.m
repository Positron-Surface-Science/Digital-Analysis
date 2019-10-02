function [filtered,d,baseline] = bpfilter(input,c,FFTIn,d)

% --------------------------------------------------------------------------
% Fast Fourier Transform band pass filter.
% --------------------------------------------------------------------------
% Should be converted to use short-time Fourier transform.
% --------------------------------------------------------------------------

% FFT filter constants.

Fs = input.desc.fs; % Sampling frequency.
N = numel(input.y); % Number of samples in waveform.
dF = Fs/N; % Samples per second / number of samples = frequency.
f = (-Fs/2:dF:Fs/2-dF)'; % Full frequency spectrum.
timebase = input.info.TIMEBASE;
center = timebase*5;
offset = -input.info.OFFSET;
%BPF = ((FFTIn(c,2) < abs(f)) & (abs(f) < FFTIn(c,1)));%0.50*10.^8));

%{
average1 = mean(input.y(100:500));
%
average2 = mean(input.y(600:1000));

lastwarn('');
pFit = polyfit([average1 average2],[input.x(250) input.x(650)],1);
[warnMsg, ~] = lastwarn;

if isempty(warnMsg)
    slope = pFit(1);
    intercept = pFit(2);
    
else
    slope = 0;
    intercept = 0;
    
end
%}

% ---------------------------------
% Try the continuous wave transform for both the band-pass filter and the
% neural network input.
% ---------------------------------

try
    baseline = mean(input.y(round((offset-999E-9)*Fs):round((offset-500E-9)*Fs)));
catch
    baseline = 0;
end

%assignin('base','baseline',baseline);
%evalin('base','baselines = horzcat(baselines,baseline);');

signal = input.y - baseline; %(average1 + average2)/2 - input.x*slope;
%figure
%[spectrum,~,~] = spectrogram(input.y,[],[],[],Fs,'yaxis');
%plot(f,real(spectrum))
%figure
%[p,f] = periodogram(signal(100:3000),[],f,Fs);
%plot(f,p);

if isempty(d) || isempty(d{c})
    %{
    d{c} = designfilt('bandpassiir','HalfPowerFrequency1',FFTIn(c,2), ...
        'HalfPowerFrequency2',FFTIn(c,1), 'FilterOrder',10,...
        'SampleRate',Fs,'DesignMethod','butter');
    %}
    %{
    d{c} = designfilt('bandpassiir','PassbandFrequency1',FFTIn(c,2), ...
        'PassbandFrequency2',FFTIn(c,1), 'FilterOrder',10,...
        'SampleRate',Fs,'PassbandRipple',0.5,'DesignMethod','cheby1');
        %}
        
        d{c} = designfilt('bandpassiir','PassbandFrequency1',FFTIn(c,2), ...
            'PassbandFrequency2',FFTIn(c,1), 'FilterOrder',20, ...
            'SampleRate',Fs,'PassbandRipple',0.25,'DesignMethod','ellip', ...
        'StopbandAttenuation1',10,'StopbandAttenuation2',10);
    
    %{
   d{c} = designfilt('bandpassfir','CutoffFrequency1',FFTIn(c,2), ...
            'CutoffFrequency2',FFTIn(c,1), 'FilterOrder',10, ...
            'SampleRate',Fs);
    %}
    %sys = mkfilter(17000,1,'rc');
end

%{
[y,t] = impulse(sys);
assignin('base','y',y);
assignin('base','t',t);
evalin('base','plot(ax2,t,y);');
%}

%{
[h,t] = stepz(d{c});
[h2,t2] = impz(d{c});
[h3,w] = freqz(d{c});
assignin('base','h',h);
assignin('base','t',t);
evalin('base','plot(ax2,t,h);');
assignin('base','h2',h2);
assignin('base','t2',t2);
evalin('base','plot(ax3,t2,h2);');
assignin('base','h3',h3);
assignin('base','w',w);
evalin('base','plot(ax4,w,h3);');
%}

%[cA,cD] = dwt(signal,'db4');

%[wt,f] = cwt(signal,Fs,'FrequencyLimits',[6E3 1E8]);

%assignin('base','s',s);

%spectrum = fftshift(fft(signal))/N;

%spectrumlpf = zeros(size(spectrum));
%size(spectrumlpf)
%cwt(signal,Fs);

%spectrumlpf = BPF.*spectrum;

%filtered = bandpass(signal,[FFTIn(c,2) FFTIn(c,1)],Fs,'ImpulseResponse','iir','Steepness',0.5);

if c == 3
    % ----------------------------------------------------
    % Neural Network stuff
    % ----------------------------------------------------
    
    %signal = input.y - mean(input.y(100:800)); %(average1 + average2)/2 - input.x*slope;
    %spectrum = fft(signal))/N;
    %spectrum = ifftshift(spectrum);
    
    %{
    figure
    plot(imag(spectrum(1:200)));
    filtered = [];
    %}
    
    filteredR = real(spectrumlpf);
    filteredI = imag(spectrumlpf);
    indexR = (numel(filteredR)/2) + 1;
    indexI = (numel(filteredI)/2) + 1;
    
    filteredR = filteredR(indexR:indexR+499);
    filteredI = filteredI(indexI:indexI+499);
    
    %{
    figure
    plot(filteredR)
    figure
    plot(filteredI)
    %}
    
    filtered = vertcat(filteredR,filteredI);
    %figure
    %plot(filteredR);
    %assignin('base','filtered',filtered);
    %evalin('base','dataVectorFFTDL = horzcat(dataVectorFFTDL,filtered);');
%}
%assignin('base','spectrumlpf',real(spectrumlpf(30500:32000)));
%evalin('base','dataVectorFFT2 = horzcat(dataVectorFFT2,spectrumlpf);');
%filtered = ifft(ifftshift(spectrumlpf))*N;
%filtered = filtered - mean(filtered(100:800));

%maximum = max(filtered);
%assignin('base','maximum',maximum);
%evalin('base','dataVectorPreshaped2 = horzcat(dataVectorPreshaped2,maximum);');
%filtered = ifft(ifftshift(spectrumlpf))*N;
%filtered = filtered - mean(filtered(100:800));

else
    
    %filtered = wdenoise(input.y,10);
    
    %filtered = ifft(ifftshift(spectrumlpf))*N;
    
    %filtered = filtered - mean(filtered(200:800));
    
    filtered = filter(d{c},signal);
    %filtered = idwt(cA,cD,'db4');
    
    %filtered = icwt(wt,f,[FFTIn(c,2) FFTIn(c,1)]);
    
end

end