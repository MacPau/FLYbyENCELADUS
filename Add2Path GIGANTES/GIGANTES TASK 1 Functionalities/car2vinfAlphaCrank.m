function [NODE] = car2vinfAlphaCrank(vinfCAR, SVMoon)


rfb= SVMoon(1:3);
vfb= SVMoon(4:6);
norm_vfb   = norm(vfb);

if isequaln( vinfCAR, zeros(1,3) )
    normVinf = NaN;
    pump     = NaN;
    crank    = NaN;
    return
end

tangential  = vfb / norm_vfb;
cross_track = cross(rfb, vfb);
cross_track = cross_track / norm(cross_track);
normal      = cross(tangential, cross_track);

vinf_tcn_cart = [dot(vinfCAR, tangential), dot(vinfCAR, cross_track), dot(vinfCAR, normal)];

normVinf = norm(vinf_tcn_cart);
pump     = acos(vinf_tcn_cart(1) / normVinf);

if abs(pump) >= 0.001
    crank = -atan2(vinf_tcn_cart(2), vinf_tcn_cart(3));
else
    crank = 0.0;
end

NODE=[normVinf, pump, crank]; 

end