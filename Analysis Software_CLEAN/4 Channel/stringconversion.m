function s = stringconversion(n)
s = [];
   while n > 0
      d = mod(n,10);
      s = [char(48+d), s];
      n = (n-d)/10;
   end
end