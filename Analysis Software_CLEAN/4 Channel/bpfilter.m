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

% ---------------------------------
% Try the continuous wave transform for both the band-pass filter and the
% neural network input.
% ---------------------------------

try
    baseline = mean(input.y(round((offset-999E-9)*Fs):round((offset-500E-9)*Fs)));
catch
    baseline = 0;
end

signal = input.y - baseline;

if isempty(d) || isempty(d{c})
    
        d{c} = designfilt('bandpassiir','PassbandFrequency1',FFTIn(c,2), ...
            'PassbandFrequency2',FFTIn(c,1), 'FilterOrder',20, ...
            'SampleRate',Fs,'PassbandRipple',0.25,'DesignMethod','ellip', ...
        'StopbandAttenuation1',10,'StopbandAttenuation2',10);
    
end

if c == 3
    
    filteredR = real(spectrumlpf);
    filteredI = imag(spectrumlpf);
    indexR = (numel(filteredR)/2) + 1;
    indexI = (numel(filteredI)/2) + 1;
    
    filteredR = filteredR(indexR:indexR+499);
    filteredI = filteredI(indexI:indexI+499);
    filtered = vertcat(filteredR,filteredI);
    
else
    
    filtered = filter(d{c},signal);
    
end

end