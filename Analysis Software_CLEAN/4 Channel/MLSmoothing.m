function [smoothPulse,error,a] = MLSmoothing(pulse1,Ts)

pulse = pulse1(round(10E-6/Ts)+240:round(30E-6/Ts)-1084)';
assignin('base','vIn',pulse);
pulse2 = pulse;


s = evalin('base','s;');

%s = s(13E-6/Ts:22E-6/Ts);

sum1 = 0;
sum2 = 0;
a = 0;
b = 0;

for it=1:5
    for n=1:numel(pulse-1)
        
        sum1 = sum1 + pulse(n)*s(n) - b*s(n);
        sum2 = sum2 + s(n)^2;
        
    end
    
    a = sum1/sum2;
    
    sum1 = 0;
    sum2 = 0;
    
    for n=1:numel(pulse-1)
        
        sum1 = sum1 + pulse(n) - a*s(n);
        
    end
    
    b = sum1/numel(pulse);
    
    sum1 = 0;
    sum2 = 0;
    
end

sum1 = 0;

x = a*s + b;

a = a;% + b;
%baseline = mean(x(250:300))

error = immse(x,pulse);
smoothPulse = x';%vertcat(pulse1,x',pulse1);
%plot([smoothPulse pulse1])
end