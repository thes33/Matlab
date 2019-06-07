function pc = invRau(r)
% Inverse RAU - rationalized arcsine transform
%

if any(r < -23 | r > 123)
   error('r should be between -23 and 123');
end

t = (r+23)/46.47324337;
pc = (sin(t/2)).^2;