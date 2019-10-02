Version differences between SH and CPB:
1)	PEGS_4channel.mlapp - line  308 in code view-functions(multiStopping)
	if isempty(timeOfFlight2{s}) == 0 && (timeOfFlight2{s} > 6000E-9 || ...
	
2)	paes.m code - Integration Filter line 119-136:
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
	
3)	tripleCoincidence_GeTiming_Trapezoidal_Functional_Input.m -	Oscilloscope Filenameing convention:
	line 54: if (true)
	line 71:  if (false)

PEGS-4Channel.mlapp Settings:
Channel Configuration ch 1: (MCP) Bining-2126
Channel Configuration ch 2: (Gamma) ELET and moving average semi gaussian shaping type (4096 bins)
Band-Pass Filter: ON; ch1 (MCP) 4-200 MHz; ch2 (Gamma) 0.0001 to 10 MHz 
Shaping/Fitting ch2: (Gamma) 5 us(top) and 2 us(rise/fall)
ToF Spectrum Controls ch1: (MCP) 4 us timing window, 
Gamma Spectrum Controls ch2: (Gamma) gamma window 12



Multi-stop Conditions Explanations:
1 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY 1 MCP pulse
2 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY 2 MCP pulse
3 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY 3 MCP pulse
All: Constructs 1 PAES spectrum using every MCP pulse (1,2,3, etc.) in each MCP trace
> 1 Only: Constructs 1 PAES spectrum using MCP traces containing ONLY MORE THAN 1 MCP pulse (2,3, etc.); if ONLY 1 MCP pulse appears in the MCP trace it throws it out
2 Only - Split Spectra: Contructs 2 PAES spectra using MCP pulses containing ONLY 2 MCP pulses
3 Only - Split Spectra: Contructs 3 PAES spectra using MCP pulses containing ONLY 3 MCP pulses