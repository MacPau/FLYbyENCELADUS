function vinfTCN = vinfAlphaCrank_to_VinfTCN(vinf, alpha, k)

% --> from (vinf, alpha, k) to (vinfT, vinfC, vinfN)

vinfTCN = vinf.*[ cos(alpha), -sin(alpha).*sin(k), sin(alpha).*cos(k) ];

end