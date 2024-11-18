function [mu, rcm, mum, radm, hmin] = centralBody_Planet(idpl, pl)

if nargin == 2

    if idpl == 1 % --> central body is SUN
        mu = 132724487690;
        [mum, radm, rcm] = planetConstants(pl);
        [hmin] = maxmin_flybyAltitude(pl);
    elseif idpl == 5 % --> central body is JUPITER
        mu  = planetConstants(idpl);
        [rcm, mum, radm, hmin] = jupMoonsConstants(pl);
    elseif idpl == 6 % --> central body is SATURN
        mu  = planetConstants(idpl);
        [rcm, mum, radm, hmin] = satMoonsConstants(pl);
    elseif idpl == 7
        mu  = planetConstants(idpl);
        [rcm, mum, radm, hmin] = uranusMoonsConstants(pl);
    end

elseif nargin == 1
    
    if idpl == 1
        mu = 132724487690;
    elseif idpl == 5
        mu  = planetConstants(idpl);
    elseif idpl == 6
        mu  = planetConstants(idpl);
    elseif idpl == 7
        mu  = planetConstants(idpl);
    end

    rcm  = [];
    mum  = [];
    radm = [];
    hmin = [];

end

end