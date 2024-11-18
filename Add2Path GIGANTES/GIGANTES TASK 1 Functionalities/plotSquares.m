function fig = plotSquares(INPUT, holdon)

if nargin == 1
    holdon = 0;
elseif nargin == 2
    if isempty(holdon)
        holdon = 0;
    end
end

if holdon == 1
    fig = gcf;
else
    fig = figure('Color', [1 1 1]);
end

hold on;

rect = INPUT.rect;

for indr = 1:length(rect)
    plot( rect(indr).rect(:,1), rect(indr).rect(:,2), 'Color', 'black',...
        'LineWidth', 2, 'HandleVisibility', 'Off' );
end

end
