% from https://www.mathworks.com/help/images/ref/imsegkmeans.html

% for 20190520 C2 specially. Finally, dropped.

% Do gray scale segmentation of ions in WideFiled image by kmeans classification.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 8;



%%
% pattern = zeros(Frame, 1);
% % Get the name of the image the user wants to use.
% folder = pwd;
% filename = 'kmeans_grayscale_brain.png';
% % [filename, folder] = uigetfile({'*.jpg';'*.bmp'},'Select image');
% fullFileName = fullfile(folder, filename);
% grayImage = imread(fullFileName);
% for ii = 1:Frame
ii = 300;
grayImage = tif{ii, 1};
grayImage = grayImage.*(~Mask1);
grayImage = imboxfilt(grayImage, 7); % Filter is important here.
% Normalize to uint8
ymax = 255; ymin=0;
xmax = max(max(grayImage)); %求得InImg中的最大值
xmin = min(min(grayImage)); %求得InImg中的最小值
grayImage = round((ymax-ymin)*(grayImage-xmin)/(xmax-xmin) + ymin);



% Display the image
subplot(2, 2, 1);
imshow(grayImage, 'DisplayRange',[], 'InitialMagnification', 'fit');
% axis on image;
% [folder, baseFileNameNoExt, ext] = fileparts(fullFileName);  % Do again in case user changes/alters the code to use a different filename.
% caption = sprintf('Input Image: "%s"', baseFileNameNoExt);
% title(caption, 'fontsize', fontSize, 'Interpreter', 'none');

% Convert to gray scale if needed.
[rows, columns, numberOfColorChannels] = size(grayImage);
if numberOfColorChannels == 3
	fprintf('That was a color image.  I am converting it to grayscale.\n');
	grayImage = rgb2gray(grayImage);
end

% Set up figure properties:
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% % Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
% 
% drawnow;
% hp = impixelinfo(); % Set up status line to see values when you mouse over the image.

% Define some number of clusters that you know will definitely be there.
numberOfClusters = 3;

% Do kmeans clustering on the gray scale image.
grayLevels = double(grayImage(:)); % Convert to column vector.
[clusterIndexes, clusterCenters] = kmeans(grayLevels, numberOfClusters,...
	'distance', 'sqEuclidean', ...
	'Replicates', 2);
labeledImage = reshape(clusterIndexes, rows, columns);
subplot(2, 2, 2);
imshow(labeledImage, 'DisplayRange',[], 'InitialMagnification', 'fit')
% caption = sprintf('k-means with %d clusters', numberOfClusters);
% title(caption,'FontSize', fontSize);
% axis on image;
% colorbar
% hp = impixelinfo(); % Set up status line to see values when you mouse over the image.

% Optional - for fun and to distinguish the classes better by giving each a unique color.
% Let's assign each blob a different color to visually show the user the distinct blobs.
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
% coloredLabels is an RGB image.  We could have applied a colormap instead (but only with R2014b and later)
subplot(2, 2, 3);
imshow(coloredLabels, 'DisplayRange',[], 'InitialMagnification', 'fit');
% axis on image; % Make sure image is not artificially stretched because of screen's aspect ratio.
% title('Pseudo colored labels, from label2rgb().', 'FontSize', fontSize);

%===============================================================================
% Now kmeans can give a different index to the tumor in each run,
% so we'll assume the tumor is the brightest class.
% Find the brightest class.  The cluster center will be the mean gray level of the class.
[minValue, indexOfMinValue] = min(clusterCenters);
% fprintf('The tumor is class #%d, with a mean gray level of %.2f\n', indexOfMaxValue, minValue);
% Get pixels that are labeled as the tumor.
tumor = labeledImage == indexOfMinValue;
% Extract the largest blob;
tumor = bwareafilt(tumor, 1);
% Fill holes.
tumor = imfill(tumor, 'holes');
% Display the image.
subplot(2, 2, 4);
imshow(tumor, 'DisplayRange',[], 'InitialMagnification', 'fit');
% axis on image;
% title('Tumor Binary Mask Image','FontSize', fontSize);
% hp = impixelinfo(); % Set up status line to see values when you mouse over the image.

% pattern(ii, 1) = sum(tumor(:));
% s = sum(tumor(:));
% if ii == 1 
% pattern(ii, 1) = s;
% else
%     if s < pattern(ii-1, 1)
%         pattern(ii, 1) = pattern(ii-1, 1);
%     else
%         pattern(ii, 1) = s;
%     end
% end
% 
% end
%%
D = (pattern(2:end, 1) - pattern(1:end-1, 1))*125*125*10^(-18)/0.01;
time = (917+1:917+Frame-1)'/100; % from 
figure('color', 'w'); 
plot(time, D)