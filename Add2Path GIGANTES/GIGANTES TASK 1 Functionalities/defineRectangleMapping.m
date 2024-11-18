function rect = defineRectangleMapping()

% --> very high latitudes [1,4] --> score: 2
rect(1).rect      = [ -180 70; -90 70; -90 90; -180 90; -180 70 ];
rect(2).rect      = rect(1).rect;
rect(2).rect(:,1) = rect(2).rect(:,1) + 90;
rect(3).rect      = rect(2).rect;
rect(3).rect(:,1) = rect(3).rect(:,1) + 90;
rect(4).rect      = rect(3).rect;
rect(4).rect(:,1) = rect(4).rect(:,1) + 90;

% --> mid-high latitudes [20,23] --> score: 1
rect(20).rect      = [ -180 30; -90 30; -90 70; -180 70; -180 30 ];
rect(21).rect      = rect(20).rect;
rect(21).rect(:,1) = rect(21).rect(:,1) + 90;

rect(22).rect      = rect(21).rect;
rect(22).rect(:,1) = rect(22).rect(:,1) + 90;

rect(23).rect      = rect(22).rect;
rect(23).rect(:,1) = rect(23).rect(:,1) + 90;

% --> central latitudes [5,10&15] --> score: 2
rect(5).rect      = [-150 -30; -90 -30; -90 30; -150 30; -150 -30];
rect(6).rect      = rect(5).rect;
rect(6).rect(:,1) = rect(6).rect(:,1) + 60;
rect(7).rect      = rect(6).rect;
rect(7).rect(:,1) = rect(7).rect(:,1) + 60;
rect(8).rect      = rect(7).rect;
rect(8).rect(:,1) = rect(8).rect(:,1) + 60;
rect(9).rect      = rect(8).rect;
rect(9).rect(:,1) = rect(9).rect(:,1) + 60;

rect(10).rect = [ 150 -30; 180 -30; 180 30; 150 30; 150 -30 ];
rect(15).rect = [ -180 -30; -150 -30; -150 30; -180 30; -180 -30 ]; % --> consider it as the same of face 10

% --> mid-low latitudes [11,14] --> score: 3
rect(11).rect      = [-180 -70; -90 -70; -90 -30; -180 -30; -180 -70 ];
rect(12).rect      = rect(11).rect;
rect(12).rect(:,1) = rect(12).rect(:,1) + 90;

rect(13).rect      = rect(12).rect;
rect(13).rect(:,1) = rect(13).rect(:,1) + 90;
rect(14).rect      = rect(13).rect;
rect(14).rect(:,1) = rect(14).rect(:,1) + 90;

% --> very-low latitudes [16,19] --> score: 4
rect(16).rect      = [-180 -90; -90 -90; -90 -70; -180 -70; -180 -90 ];
rect(17).rect      = rect(16).rect;
rect(17).rect(:,1) = rect(17).rect(:,1) + 90;

rect(18).rect      = rect(17).rect;
rect(18).rect(:,1) = rect(18).rect(:,1) + 90;
rect(19).rect      = rect(18).rect;
rect(19).rect(:,1) = rect(19).rect(:,1) + 90;

for indr = 1:length(rect)
    rect(indr).face = indr;
    if indr >= 1 && indr <= 4
        rect(indr).faceValue = 2;
    elseif indr >= 20 && indr <= 23
        rect(indr).faceValue = 1; % --> very high value for low-latitude
    elseif indr >= 5 && indr <= 10 % --> high value for low-latitude
        rect(indr).faceValue = 2;
    elseif indr >= 11 && indr <= 14 % --> mid value for mid-latitude
        rect(indr).faceValue = 3;
    elseif indr >= 16 && indr <= 19 % --> low value for high latitude
        rect(indr).faceValue = 4;
    elseif indr == 15
        rect(indr).faceValue = 2;
    end
end

% --> make it comparable with Campagnola --> shift to east longitude
for indr = 1:length(rect)
    rect(indr).rect(:,1) = rect(indr).rect(:,1) + 180;
end

end
