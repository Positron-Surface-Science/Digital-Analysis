paes.m code:
Integration Filter:
	% Integrating peak and region before peak to reduce periodic noise
		% pulses.
		if TypeIn == 1 && VMinIndex{c}(s)-round((10*10.^(-9)/dN)) > 0
			integration = trapz(pulseIny(VMinIndex{c}(s)-round((10*10.^(-9)/dN)) ...
				:VMinIndex{c}(s)));
		elseif TypeIn == 1
			integration = 0;
			%'TypeIn 1'
		end
		
		if TypeIn == 1 && integration < 0.025
		   %numPeaks = numPeaks - 1;
		   crossTime(s) = NaN;
		   'integration'
		   integration
		   continue;
			
		end
Polarity of MCP pulse and peak finding function:
	if TypeIn == 1
    pulseIny = -pulseIny;
    lastwarn('');
    % Set to detect peaks 15 ns apart (15*10.^(-9)/dN)
  [VMin{c},VMinIndex{c}] = findpeaks(pulseIny(50:numel(pulseIny)-50),...
        'MinPeakDistance',round((15*10.^(-9))/dN), ...
        'MinPeakHeight',2*10.^(-3),'MinPeakProminence',2*10.^(-3), ...
        'MinPeakWidth',round((2*10.^(-9))/dN),'WidthReference', ...
        'halfheight','MaxPeakWidth',round((7*10.^(-9)/dN)));
		
tripleCoincidence_GeTiming_Trapezoidal_Functional_Input.m:
Oscilloscope Filenameing convention:
	line 54: true
	line 71: false

PEGS-4Channel.mlapp Settings:
Bining-2126/4096
band-pass filter-ch1: 4-200 MHz; ch2: 0.0001 to 10 MHz 
shaping/fitting-ch2: 5 us(top) and 2 us(rise/fall)
ToF Spectrum Controls: 4 us timing window, 
Gamma Spectrum Contros: ch2: gamma window 12
Channel Configuration: ch 2 (Gamma) is ELET and moving average semi gaussian shaping type


Multi-stop Conditions Explanations:
1 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY 1 MCP pulse
2 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY 2 MCP pulse
3 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY 3 MCP pulse
All: Constructs 1 PAES spectrum using every MCP pulse (1,2,3, etc.) in each MCP trace
> 1 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY MORE THAN 1 MCP pulse (2,3, etc.); if ONLY 1 MCP pulse appears in the MCP trace it throws it out
2 Only - Split Spectra: Contructs 2 PAES spectra using MCP pulses containing ONLY 2 MCP pulses
3 Only - Split Spectra: Contructs 3 PAES spectra using MCP pulses containing ONLY 3 MCP pulses