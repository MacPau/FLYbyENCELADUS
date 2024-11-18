function [alpha, vinf] = raRp2AlphaVinf_v2(ra, rp, ADIM, mu)

% on a Tisserand plot, this function computes (alpha,vinf) from (ra,rp)
% only valid for elliptical orbits

% moon properties
rPL = ADIM(1);
vPL = sqrt(mu/rPL);

a = 0.5.*(ra + rp)./rPL; % non-dimensional semi-major axis
e = (ra - rp)./(ra + rp);

% compute the corresponding vinf
vinf = sqrt(3 - 1./a - sqrt(4.*a.*(1 - e.^2))).*vPL; % dimensional vinf (km/s)

% compute the corresponding alpha
alpha = acos(vPL./(2.*vinf).*(1 - (vinf./vPL).^2 - 1./a));
alpha = wrapToPi(alpha);

end