function [DIFF, DV, tofsc, node, alpha1, asc1, asc2, esc1, esc2, rla, th1, th2, T1, T2] =...
    solveForAL1(alpha1, ADIM, S, vinf1, vinf2, mu)

% DESCRIPTION : 
%
% this function is used to compute VILTs for a given Saturnian moon, given
% intial and arrival orbit on a Tisserand map (i.e. vinf1 and vinf2 need to
% be specified). A recursive procedure is implemented to find alpha1.
% 
% INPUT : 
% 
% alpha1 : post-flyby pump angle (unknown) --> a recursive procedure is
%          used to find it
% ADIM   : 1x3 vector where ADIM(1)=rcm (circular orbit radius of the moon)
% S      : 1x5 vector containing variables defining the VILT
%          S(1) = TYPE --> OO/OI/IO/II (88, 81, 18, 11)
%          S(2) = kei  --> +1 for EXTERNAL, -1 for INTERNAL
%          S(3) = N    --> number of moon revolutions
%          S(4) = M    --> number of SC revolutions
%          S(5) = L    --> number of SC revolution where the DSM occurs
% vinf1  : initial infinity velocity w.r.t. the moon
% vinf2  : final infinity velocity w.r.t. the moon
%
% OUTPUT : 
% 
% DIFF   : difference between SC and moon transfer times (to solve the
%          phasing problem)
% DV     : non-dimensional DV of the VILT
% tofsc  : non-dimensional TOF of the SC on the VILT
% node   : 1x2 vector containing the new node (alpha,vinf)
% alpha1 : post-flyby pump angle (rad)
% asc1   : non-dimensional SC semi-major axis of the orbit before the manoeuvre
% asc2   : non-dimensional SC semi-major axis of the orbit after the manoeuvre
%          (arrival node)
% esc1   : SC eccentricity of the orbit before the manoeuvre
% esc2   : SC eccentricity of the orbit after the manoeuvre (arrival node)
% rla    : non-dimensional leveraging apse (i.e. where the manoeuvre is located,
%          apoapsis or periapsis?)
% th1    : true anomaly at the first moon encounter (rad)
% th2    : true anomaly at the second moon encounter (rad)
%
% -------------------------------------------------------------------------

if size(alpha1,1) == 1
    alpha1 = alpha1';
end

type    = S(1);
int_ext = S(2);
N       = S(3);
M       = S(4);
L       = S(5);

% solve for alpha1
[DIFF, DV, tofsc, alpha1, asc1, asc2, esc1, esc2, rla, th1, th2, T1, T2] = ...
    solveForAlpha1(alpha1, ADIM, int_ext, type, N, M, L, vinf1, vinf2, mu);

% find the new node on Tisserand map (alpha,vinf)
rp2            = asc2.*(1 - esc2)*ADIM(1); % dimensional rp
ra2            = asc2.*(1 + esc2)*ADIM(1); % dimensional ra
[alpha2,vinf2] = raRp2AlphaVinf_v2(ra2, rp2, ADIM, mu);
node           = [alpha2 vinf2];

if ~isreal(DIFF)
    for inddif = 1:length(DIFF)
        if imag(DIFF(inddif)) < 1e-6
            DIFF(inddif) = real(DIFF(inddif));
        end
    end
end

end