function [row, col] = ImageJroiLocation(sROI)
switch sROI.strType
    case 'Rectangle'
        RectBounds = sROI.vnRectBounds;
        % [x-left_up y-left_up x-right_down y-right_down]
        col = [RectBounds(2), RectBounds(4)];
        row = [RectBounds(1), RectBounds(3)];
    case 'Polygon'
        Polygon = sROI.mnCoordinates;
        col = [min(Polygon(:,2)), max(Polygon(:,2))];
        row = [min(Polygon(:,1)), max(Polygon(:,1))];
    case 'Oval'
        RectBounds = sROI.vnRectBounds;
        % [x-left_up y-left_up x-right_down y-right_down]
        col = [RectBounds(2), RectBounds(4)];
        row = [RectBounds(1), RectBounds(3)];
end