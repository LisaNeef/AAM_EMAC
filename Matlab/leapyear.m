function L = leapyear(yy)
%--------leapyear.m
% given an input array, determine which years are leapyears
% here I'm using the algorithm given on wikipedia:
%if year modulo 400 is 0
%       then is_leap_year
%else if year modulo 100 is 0
%       then not_leap_year
%else if year modulo 4 is 0
%       then is_leap_year
%else
%       not_leap_year


L = yy*0;

for iy = 1:length(yy)
  if mod(yy(iy),400) == 0
     L(iy) = 1;
  else
     if (mod(yy(iy),100) ~= 0) & (mod(yy(iy),4) == 0)
       disp(['leap year!  ' num2str(yy(iy))]);
       L(iy) = 1;
     end
  end
end


end
