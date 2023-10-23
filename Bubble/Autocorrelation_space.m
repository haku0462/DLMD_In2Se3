%%
% Read the polarization distribution image
polarizationImage = imread('11-0-250.png');

% Convert the image to binary
binaryImage = im2bw(polarizationImage, graythresh(polarizationImage));

% Perform domain segmentation using region-growing
labeledImage = bwlabel(binaryImage);
stats = regionprops(labeledImage, 'Area');

% Extract domain sizes
domainSizes = [stats.Area];

% Compute autocorrelation of the binary image
autocorrImage = normxcorr2(binaryImage, binaryImage);

% Find the peak in the autocorrelation function
[~, maxIdx] = max(autocorrImage(:));
[maxRow, maxCol] = ind2sub(size(autocorrImage), maxIdx);

% Calculate the average domain size
avgDomainSize = mean(domainSizes);

% Display the results
figure;
subplot(1, 2, 1);
imshow(polarizationImage);
title('Polarization Distribution');

subplot(1, 2, 2);
imshow(autocorrImage);
hold on;
plot(maxCol, maxRow, 'r+', 'MarkerSize', 10);
hold off;
title('Autocorrelation');
colorbar;

fprintf('Average Domain Size: %.2f pixels\n', avgDomainSize);
fprintf('Spatial Distribution Peak Location: Row %d, Column %d\n', maxRow, maxCol);
%%
% Read the polarization distribution image
polarizationImage = imread('11-0-250.png');

% Convert the image to binary
binaryImage = im2bw(polarizationImage, graythresh(polarizationImage));

% Invert the binary image to set the dots as foreground
binaryImage = ~binaryImage;

% Perform domain segmentation using connected component analysis
cc = bwconncomp(binaryImage);
numPixels = cellfun(@numel, cc.PixelIdxList);
[~, idx] = max(numPixels);
binaryImage = false(size(binaryImage));
binaryImage(cc.PixelIdxList{idx}) = true;

% Perform morphological operations to enhance the domain regions
se = strel('disk', 3);
binaryImage = imopen(binaryImage, se);

% Perform domain segmentation again to obtain separate domains
labeledImage = bwlabel(binaryImage);
stats = regionprops(labeledImage, 'Area');

% Extract domain sizes
domainSizes = [stats.Area];

% Compute autocorrelation of the binary image
autocorrImage = normxcorr2(double(binaryImage), double(binaryImage));

% Find the peak in the autocorrelation function
[~, maxIdx] = max(autocorrImage(:));
[maxRow, maxCol] = ind2sub(size(autocorrImage), maxIdx);

% Calculate the percentage of domain pixels
domainPixels = sum(binaryImage(:));
totalPixels = numel(binaryImage);
domainPercentage = (domainPixels / totalPixels) * 100;

% Display the results
figure;
subplot(2, 2, 1);
imshow(polarizationImage);
title('Polarization Distribution');

subplot(2, 2, 2);
imshow(binaryImage);
title('Polarization Domain');

subplot(2, 2, 3);
imshow(autocorrImage);
hold on;
plot(maxCol, maxRow, 'r+', 'MarkerSize', 10);
hold off;
title('Autocorrelation');
colorbar;

fprintf('Average Domain Size: %.2f pixels\n', mean(domainSizes));
fprintf('Domain Percentage: %.2f%%\n', domainPercentage);
%%
% Read the polarization distribution image
polarizationImage = imread('11-40-800.png');

% Convert the image to grayscale
grayImage = rgb2gray(polarizationImage);

% Threshold the image to separate the red and blue dots
redThreshold = 50; % Adjust this threshold value for red dots
blueThreshold = 100; % Adjust this threshold value for blue dots
binaryRedDots = grayImage > redThreshold;
binaryBlueDots = grayImage > blueThreshold;

% Perform morphological operations to enhance the dots
se = strel('disk', 2);
binaryRedDots = imclose(binaryRedDots, se);
binaryBlueDots = imclose(binaryBlueDots, se);

% Perform connected component analysis to identify the red dot fields
ccRed = bwconncomp(binaryRedDots);
numPixelsRed = cellfun(@numel, ccRed.PixelIdxList);

% Adjust thresholdSize for field size
thresholdSize = 100; % Adjust this value for field size
fieldIndicesRed = find(numPixelsRed > thresholdSize);
binaryRedFields = false(size(binaryRedDots));
for i = 1:numel(fieldIndicesRed)
    binaryRedFields(ccRed.PixelIdxList{fieldIndicesRed(i)}) = true;
end

% Perform connected component analysis to identify the blue dot fields
ccBlue = bwconncomp(binaryBlueDots);
numPixelsBlue = cellfun(@numel, ccBlue.PixelIdxList);

% Adjust thresholdSize for field size
thresholdSize = 100; % Adjust this value for field size
fieldIndicesBlue = find(numPixelsBlue > thresholdSize);
binaryBlueFields = false(size(binaryBlueDots));
for i = 1:numel(fieldIndicesBlue)
    binaryBlueFields(ccBlue.PixelIdxList{fieldIndicesBlue(i)}) = true;
end

% Calculate the sizes of red and blue fields
redFieldSizes = numPixelsRed(fieldIndicesRed);
blueFieldSizes = numPixelsBlue(fieldIndicesBlue);

% Compute the percentage of red and blue fields
totalPixels = numel(polarizationImage);
redFieldPercentage = sum(redFieldSizes) / totalPixels * 100;
blueFieldPercentage = sum(blueFieldSizes) / totalPixels * 100;

% Compute the autocorrelation of the red and blue fields
autocorrRed = normxcorr2(double(binaryRedFields), double(binaryRedFields));
autocorrBlue = normxcorr2(double(binaryBlueFields), double(binaryBlueFields));

% Find the peak in the autocorrelation functions
[~, maxIdxRed] = max(autocorrRed(:));
[maxRowRed, maxColRed] = ind2sub(size(autocorrRed), maxIdxRed);
[~, maxIdxBlue] = max(autocorrBlue(:));
[maxRowBlue, maxColBlue] = ind2sub(size(autocorrBlue), maxIdxBlue);

% Display the results
% figure;
% subplot(2, 2, 1);
% imshow(polarizationImage);
% title('Polarization Distribution');
% 
% subplot(2, 2, 2);
% imshow(binaryRedFields);
% title('Red Fields');
% 
% subplot(2, 2, 3);
% imshow(binaryBlueFields);
% title('Blue Fields');
% 
% subplot(2, 2, 4);
% imshow(autocorrRed);
% hold on;
% plot(maxColRed, maxRowRed, 'r+', 'MarkerSize', 10);
% hold off;
% title('Autocorrelation (Red Fields)');
% colorbar;


imshow(autocorrRed);
hold on;
plot(maxColRed, maxRowRed, 'r+', 'MarkerSize', 20);
hold off;
title('Autocorrelation (Red Fields)');
colorbar;

fprintf('Red Field Percentage: %.2f%%\n', redFieldPercentage);
fprintf('Blue Field Percentage: %.2f%%\n', blueFieldPercentage);
%%
% Read the polarization distribution image
polarizationImage = imread('15-60-750.png');

% Convert the image to grayscale
grayImage = rgb2gray(polarizationImage);

% Threshold the image to separate the red and blue dots
redThreshold = 50; % Adjust this threshold value for red dots
blueThreshold = 100; % Adjust this threshold value for blue dots
binaryRedDots = grayImage > redThreshold;
binaryBlueDots = grayImage > blueThreshold;

% Perform morphological operations to enhance the dots
se = strel('disk', 2);
binaryRedDots = imclose(binaryRedDots, se);
binaryBlueDots = imclose(binaryBlueDots, se);

% Perform connected component analysis to identify the red dot fields
ccRed = bwconncomp(binaryRedDots);
numPixelsRed = cellfun(@numel, ccRed.PixelIdxList);

% Adjust thresholdSize for field size
thresholdSize = 100; % Adjust this value for field size
fieldIndicesRed = find(numPixelsRed > thresholdSize);
binaryRedFields = false(size(binaryRedDots));
for i = 1:numel(fieldIndicesRed)
    binaryRedFields(ccRed.PixelIdxList{fieldIndicesRed(i)}) = true;
end

% Perform connected component analysis to identify the blue dot fields
ccBlue = bwconncomp(binaryBlueDots);
numPixelsBlue = cellfun(@numel, ccBlue.PixelIdxList);

% Adjust thresholdSize for field size
thresholdSize = 100; % Adjust this value for field size
fieldIndicesBlue = find(numPixelsBlue > thresholdSize);
binaryBlueFields = false(size(binaryBlueDots));
for i = 1:numel(fieldIndicesBlue)
    binaryBlueFields(ccBlue.PixelIdxList{fieldIndicesBlue(i)}) = true;
end

% Calculate the sizes of red and blue fields
redFieldSizes = numPixelsRed(fieldIndicesRed);
blueFieldSizes = numPixelsBlue(fieldIndicesBlue);

% Compute the autocorrelation of the red and blue fields
autocorrRed = normxcorr2(double(binaryRedFields), double(binaryRedFields));
autocorrBlue = normxcorr2(double(binaryBlueFields), double(binaryBlueFields));

% Find the peak in the autocorrelation functions
[~, maxIdxRed] = max(autocorrRed(:));
[maxRowRed, maxColRed] = ind2sub(size(autocorrRed), maxIdxRed);
[~, maxIdxBlue] = max(autocorrBlue(:));
[maxRowBlue, maxColBlue] = ind2sub(size(autocorrBlue), maxIdxBlue);

% Enhance the visualization of autocorrelation
autocorrRed = imadjust(autocorrRed, stretchlim(autocorrRed), []);
autocorrBlue = imadjust(autocorrBlue, stretchlim(autocorrBlue), []);

% Display the results
figure;
subplot(2, 2, 1);
imshow(polarizationImage);
title('Polarization Distribution');

subplot(2, 2, 2);
imshow(binaryRedFields);
title('Red Fields');

subplot(2, 2, 3);
imshow(binaryBlueFields);
title('Blue Fields');

subplot(2, 2, 4);
imshow(autocorrRed);
hold on;
plot(maxColRed, maxRowRed, 'r+', 'MarkerSize', 10);
hold off;
title('Autocorrelation (Red Fields)');
colormap(jet); % Adjust the colormap for better visibility
colorbar;

% Save the figure as a 600 dpi image
% resolution = 600; % Adjust the resolution as desired
% filename = 'output_image';
% set(gcf, 'Units', 'inches');
% set(gcf, 'PaperPosition', [0 0 size(polarizationImage, 2) / resolution size(polarizationImage, 1) / resolution]);
% saveas(gcf, filename, 'png', '-r600');

fprintf('Red Field Percentage: %.2f%%\n', sum(redFieldSizes) / numel(polarizationImage) * 100);
fprintf('Blue Field Percentage: %.2f%%\n', sum(blueFieldSizes) / numel(polarizationImage) * 100);
%%
% Load the image
image = imread('11-40-450.png'); % Replace 'image_path.jpg' with the actual path to your image

% Convert the image to grayscale
gray = rgb2gray(image);

% Define the color thresholds for red and blue
red_threshold = 50; % Adjust the threshold value as needed
blue_threshold = 100; % Adjust the threshold value as needed

% Create masks for red and blue dots
red_mask = gray > red_threshold;
blue_mask = gray < blue_threshold;

% Apply the masks to extract red and blue dot regions
red_dots = image .* uint8(cat(3, red_mask, red_mask, red_mask));
blue_dots = image .* uint8(cat(3, blue_mask, blue_mask, blue_mask));

% Calculate the areas of red and blue dot regions
red_area = nnz(red_mask);
blue_area = nnz(blue_mask);

% Compute the autocorrelation of the red and blue dot masks
red_autocorr = normxcorr2(red_mask, red_mask);
blue_autocorr = normxcorr2(blue_mask, blue_mask);

% Set the desired value range for the colorbar
caxis_range = [-0.2, 0.5];

% Adjust the contrast of the autocorrelation images
red_min = min(red_autocorr(:));
red_max = max(red_autocorr(:));
red_autocorr_contrast = (red_autocorr - red_min) / (red_max - red_min);

blue_min = min(blue_autocorr(:));
blue_max = max(blue_autocorr(:));
blue_autocorr_contrast = (blue_autocorr - blue_min) / (blue_max - blue_min);

Plot the autocorrelation images with adjusted color mapping
subplot(1, 2, 1);
imagesc(red_autocorr_contrast);
colormap(jet);
colorbar;
title('Red Autocorrelation');

subplot(1, 2, 2);
imagesc(blue_autocorr_contrast);
colormap(jet);
colorbar;
title('Blue Autocorrelation');

% Print the calculated areas
disp(['Red Area: ', num2str(red_area)]);
disp(['Blue Area: ', num2str(blue_area)]);
% Plot the autocorrelation images
subplot(1, 2, 1);
imagesc(red_autocorr);
colormap jet;
colorbar;
title('Red Autocorrelation');

subplot(1, 2, 2);
imagesc(blue_autocorr);
colormap jet;
colorbar;
title('Blue Autocorrelation');

% Print the calculated areas
disp(['Red Area: ', num2str(red_area)]);
disp(['Blue Area: ', num2str(blue_area)]);
%%
% Load the image
image = imread('15-40-450.png'); % Replace 'image_path.jpg' with the actual path to your image

% Convert the image to grayscale
gray = rgb2gray(image);

% Define the color thresholds for red and blue
red_threshold = 50; % Adjust the threshold value as needed
blue_threshold = 100; % Adjust the threshold value as needed

% Create masks for red and blue dots
red_mask = gray > red_threshold;
blue_mask = gray < blue_threshold;

% Apply the masks to extract red and blue dot regions
red_dots = image .* uint8(cat(3, red_mask, red_mask, red_mask));
blue_dots = image .* uint8(cat(3, blue_mask, blue_mask, blue_mask));

% Calculate the areas of red and blue dot regions
red_area = nnz(red_mask);
blue_area = nnz(blue_mask);

% Compute the autocorrelation of the red and blue dot masks
red_autocorr = normxcorr2(red_mask, red_mask);
blue_autocorr = normxcorr2(blue_mask, blue_mask);

% Adjust the contrast of the autocorrelation images
red_min = min(red_autocorr(:));
red_max = max(red_autocorr(:));
red_autocorr_contrast = (red_autocorr - red_min) / (red_max - red_min);

blue_min = min(blue_autocorr(:));
blue_max = max(blue_autocorr(:));
blue_autocorr_contrast = (blue_autocorr - blue_min) / (blue_max - blue_min);

% Set the desired value range for the colorbar
caxis_range = [0.2, 0.4];

% Plot the autocorrelation images with adjusted color mapping and value range

figure('Position', [100, 100, 600, 600]);

subplot(1, 1, 1);
imagesc(red_autocorr_contrast, caxis_range);
colormap(jet);
colorbar;
title('Red Autocorrelation');



% subplot(1, 2, 2);
% imagesc(blue_autocorr_contrast, caxis_range);
% colormap(jet);
% colorbar;
% title('Blue Autocorrelation');

% Print the calculated areas
disp(['Red Area: ', num2str(red_area)]);
disp(['Blue Area: ', num2str(blue_area)]);








