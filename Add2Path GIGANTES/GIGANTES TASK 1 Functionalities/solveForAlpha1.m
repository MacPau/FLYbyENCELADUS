function [DIFF, DV, tofsc, alpha1, asc1, asc2, esc1, esc2, rla, th1, th2, T1, T2] = ...
            solveForAlpha1(alpha1, ADIM, int_ext, type, N, M, L, vinf1, vinf2, mu)

% alpha 1 is the post-flyby pump angle

rcm    = ADIM(1);

Vdim   = sqrt(mu/rcm);
vinf1  = vinf1/Vdim;
vinf2  = vinf2/Vdim;

kei    = int_ext;
kio    = type2kio(type);

asc1  = NaN.*zeros(length(alpha1),1);
rla   = NaN.*zeros(length(alpha1),1);
asc2  = NaN.*zeros(length(alpha1),1);
Tsc1  = NaN.*zeros(length(alpha1),1);
Tsc2  = NaN.*zeros(length(alpha1),1);
esc1  = NaN.*zeros(length(alpha1),1);
esc2  = NaN.*zeros(length(alpha1),1);
Esc1  = NaN.*zeros(length(alpha1),1);
Esc2  = NaN.*zeros(length(alpha1),1);
DIFF  = NaN.*zeros(length(alpha1),1);
tofsc = NaN.*zeros(length(alpha1),1);
tofga = NaN.*zeros(length(alpha1),1);
DV    = NaN.*zeros(length(alpha1),1);
th1   = NaN.*zeros(length(alpha1),1);
th2   = NaN.*zeros(length(alpha1),1);
T1    = NaN.*zeros(length(alpha1),1);
T2    = NaN.*zeros(length(alpha1),1);

for indi = 1:length(alpha1)
    
    al1          = alpha1(indi);
    asc1(indi,1) = 1/(1 - vinf1^2 - 2*vinf1*cos(al1));
    if asc1(indi,1) < 0
        
    else
        rla(indi,1)  = asc1(indi,1) + kei*sqrt(asc1(indi,1)^2 - 1/4*asc1(indi,1)*(3 - vinf1^2)^2 + 1/2*(3 - vinf1^2) - 1/(4*asc1(indi,1)));

        a = -8*rla(indi,1) + (3 - vinf2^2)^2;
        b = 4*rla(indi,1)^2 - 2*(3 - vinf2^2);
        c = 1;

        a2(1) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        a2(2) = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        a2(a2 < 0) = [];
        if ~isempty(a2)
            
            if length(a2) == 2 % ci sono due soluzioni
                a2max = max(a2);
                a2min = min(a2);
                if int_ext == +1 && ~isempty(a2((a2 < rla(indi,1)))) % exterior
                    if length(a2((a2 < rla(indi,1)))) == 1
                        asc2(indi,1) = a2((a2 < rla(indi,1)));
                    elseif length(a2((a2 < rla(indi,1)))) == 2
                        asc2(indi,1) = a2min;
                    end
                elseif int_ext == -1 && ~isempty(a2((a2 > rla(indi,1)))) % interior
                    if length(a2((a2 > rla(indi,1)))) == 1
                        asc2(indi,1) = a2((a2 > rla(indi,1)));
                    elseif length(a2((a2 > rla(indi,1)))) == 2
                        asc2(indi,1) = a2max;
                    end
                end
            else % c'è una soluzione
                if int_ext == +1 && a2 < rla(indi,1)
                    asc2(indi,1) = a2;
                elseif int_ext == -1 && a2 > rla(indi,1)
                    asc2(indi,1) = a2;
                end
            end

            if isreal(asc2(indi,1)) && ~isnan(asc2(indi,1))
                DV(indi,1) = abs(sqrt(2/rla(indi,1) - 1/asc1(indi,1)) - sqrt(2/rla(indi,1) - 1/asc2(indi,1)));

                Tsc1(indi,1) = asc1(indi,1)^(3/2); % --> period of sc orbit after the fly-by (pre-manoeuvre)
                Tsc2(indi,1) = asc2(indi,1)^(3/2); % --> period of sc orbit after the fly-by (post-manoeuvre)

                esc1(indi,1) = (1/kei)*(rla(indi,1)/asc1(indi,1) - 1); % --> eccentricity of sc orbit after the fly-by (pre-manoeuvre)
                esc2(indi,1) = (1/kei)*(rla(indi,1)/asc2(indi,1) - 1); % --> eccentricity of sc orbit after the fly-by (post-manoeuvre)

                Esc1(indi,1) = kio(1)*acos((asc1(indi,1) - 1)/(esc1(indi,1)*asc1(indi,1)));

                if isreal(Esc1(indi,1))
                    Esc2(indi,1) = kio(2)*acos((asc2(indi,1) - 1)/(esc2(indi,1)*asc2(indi,1)));

                    if isreal(Esc2(indi,1))
                        tausc1 = (Tsc1(indi,1)/(2*pi))*(Esc1(indi,1) - esc1(indi,1)*sin(Esc1(indi,1)));
                        tausc2 = (Tsc2(indi,1)/(2*pi))*(Esc2(indi,1) - esc2(indi,1)*sin(Esc2(indi,1)));

                        if type == 81
                            Ma = M + 1;
                        else
                            Ma = M;
                        end

                        tofsc(indi,1) = tausc2 - tausc1 + Tsc1(indi,1)*(L + (1 + kei)/4) ...
                            + Tsc2(indi,1)*(Ma - L - (1 + kei)/4);
                        
                        T1(indi,1) = Tsc1(indi,1)*(L + (1 + kei)/4) - tausc1;
                        T2(indi,1) = tausc2 + Tsc2(indi,1)*(Ma - L - (1 + kei)/4);

                        if type == 81
                            Na = N + 1;
                        else
                            Na = N;
                        end

                        th1(indi,1) = kio(1)*acos((1/esc1(indi,1))*(asc1(indi,1)*(1 - esc1(indi,1)^2) - 1));
                        th2(indi,1) = kio(2)*acos((1/esc2(indi,1))*(asc2(indi,1)*(1 - esc2(indi,1)^2) - 1));

                        tofga(indi,1) = Na + (1/(2*pi))*(th2(indi,1) - th1(indi,1));

                        DIFF(indi,1) = (tofga(indi,1) - tofsc(indi,1));
                    end
                end
            end
        end
    end
end

idxs         = isnan(DIFF);
DIFF(idxs)   = [];
DV(idxs)     = [];
asc1(idxs)   = [];
asc2(idxs)   = [];
esc1(idxs)   = [];
esc2(idxs)   = [];
rla(idxs)    = [];
tofsc(idxs)  = [];
alpha1(idxs) = [];
th1(idxs)    = [];
th2(idxs)    = [];
T1(idxs)     = [];
T2(idxs)     = [];

end
