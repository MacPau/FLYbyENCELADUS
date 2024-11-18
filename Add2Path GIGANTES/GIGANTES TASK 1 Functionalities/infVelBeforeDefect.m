function [vvinfouBM, delta, vvouBM] = infVelBeforeDefect(vvinfin, vvinfouAM, deltaMax, vvga)

% --> compute the infinity velocity before the DV-defect manoeuvre

delta = acos(dot(vvinfouAM, vvinfin)/(norm(vvinfin)*norm(vvinfouAM)));

if delta <= deltaMax
    vvinfouBM = norm(vvinfin).*(vvinfouAM)./norm(vvinfouAM); % --> inf. vel. out. BM
else
    b1        = vvinfin./norm(vvinfin);
    b2        = (cross(b1, vvga)./norm(vvga))./norm((cross(b1, vvga)./norm(vvga)));
    b3        = cross(b1, b2);
    MAT       = [b1' b2' b3'];

    vvInfOU   = norm(vvinfin).*(vvinfouAM)./norm(vvinfouAM);
    vec       = 1/norm(vvInfOU).*(inv(MAT)*vvInfOU');
    gam       = atan2(vec(3),vec(2)); % from Izzo
    gam       = wrapTo2Pi(gam);       % gam in [0, 360] deg
    vvinfouBM = norm(vvinfin).*[cos(deltaMax).*b1 + cos(gam)*sin(deltaMax).*b2 + sin(gam)*sin(deltaMax).*b3];
    delta = deltaMax;
end
% --> SC velocity after the flyby before the manoeuvre
vvouBM = vvinfouBM + vvga;

end