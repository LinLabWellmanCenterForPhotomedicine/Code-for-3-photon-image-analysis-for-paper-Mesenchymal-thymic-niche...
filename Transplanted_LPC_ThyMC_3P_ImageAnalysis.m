%% 

% Code for the analysis and quantification of three-photon fluorescence 
% images used in the manuscript titled: "Mesenchymal thymic niche cells 
% enable regeneration of the adult thymus and T cell immunity."
% Developed by Kamdin Mirsanaye at the Wellman Center for Photomedicine.

% Before executing the code, please ensure that the relevant image stacks 
% have been downloaded and saved from the GitHub repository.

%% Read the image stacks and remove autoflourescence

% Define the parent folder that contains the group folders
parentFolder = uigetdir;

% Check if user selected a folder or clicked cancel
if parentFolder ~= 0
    disp(['Selected folder: ', parentFolder]);
else
    disp('No folder selected.');
end

% Get a list of all .tif files in the folder
tifFiles = dir(fullfile(parentFolder, '*.tif'));  % List of .tif files in the folder

% Initialize cell arrays to hold the image stacks and folder names
imageStacks = cell(1, length(tifFiles));
folderNames_old = cell(1, length(tifFiles));

% Loop through each .tif file in the folder
for i = 1:length(tifFiles)
    
    % Get the full file path
    imageFilePath = fullfile(parentFolder, tifFiles(i).name);
    
    % Use imfinfo to get information about the TIFF file
    info = imfinfo(imageFilePath);
    
    % Check if the file has multiple pages (slices)
    numSlices = numel(info);
    
    % Initialize an array to hold the stack of images (m x n x 3 x k)
    stack = [];
    
    % Loop through each slice of the .tif file
    for slice = 1:numSlices
        % Read the current slice using imread
        currentSlice = imread(imageFilePath, slice);  % slice is the page number (1-based index)
        
        % If the current slice is RGB, it will be m x n x 3
        % Append the current slice to the stack along the 4th dimension (k)
        stack = cat(4, stack, currentSlice);  
    end
    
    % Store the stack in imageStacks
    imageStacks{i} = stack;
    
    % Extract the folder name (the parent folder of the current .tif file)
    [folderPath, folderName, ~] = fileparts(imageFilePath);
    
    % Store the folder name in the folderNames array
    folderNames_old{i} = folderName;
end

% remove outliers in control stacks
imageStacks{1} = imageStacks{1}(:,1:1000,:,:);
imageStacks{7} = imageStacks{7}(1:1000,:,:,:);

clearvars -except folderNames_old imageStacks

% Compute the channel signal ratios and Subtract the scaled channel signals

% Initialize an array to store the RGB averages for each stack
rgbAverages = zeros(3, length(imageStacks));

% Loop through each image stack
for i = 1:length(imageStacks)
    % Get the current image stack (m x n x 3 x k)
    currentStack = imageStacks{i};
    
    % Separate the RGB channels: each will have dimensions m x n x k
    R = currentStack(:,:,1,:);  % Lpc channel (m x n x k)
    G = currentStack(:,:,2,:);  % Thymc channel (m x n x k)
    B = currentStack(:,:,3,:);  % Blue channel (m x n x k)

    % Perform element-wise division between the channels
    ratioRG = double(R)./double(G);
    ratioRB = double(R)./double(B);
    ratioGB = double(G)./double(B);

    % Remove ratios that are 0 or 255 (to avoid problematic values like pure black or saturated colors)
    ratioRG(ratioRG == 0 | isinf(ratioRG) == 1) = nan;
    ratioRB(ratioRB == 0 | isinf(ratioRB) == 1) = nan;
    ratioGB(ratioGB == 0 | isinf(ratioGB) == 1) = nan;

    % Calculate the average of each ratio across all pixels (m x n) and slices (k)
    avgRG = mean(ratioRG(:),"omitnan");  % Average across all pixels and slices
    avgRB = mean(ratioRB(:),"omitnan");  % Average across all pixels and slices
    avgGB = mean(ratioGB(:),"omitnan");  % Average across all pixels and slices

    % Store the averages in the rgbAverages array (3 x n where n is the number of stacks)
    rgbAverages(:, i) = [avgRG; avgRB; avgGB];
end

% Calculate control correction factors using the first 9 stacks (assuming control samples are the first 9)
control_correction_factors = mean(rgbAverages(:, 1:8), 2);  % Mean across the first 9 columns (stacks)

% Initialize a cell array to store the cleaned image stacks
imageStacks_clean = cell(1, length(imageStacks));

% Loop through each image stack for processing
for i = 1:length(imageStacks)
    % Get the current image stack (m x n x 3 x k)
    currentStack = imageStacks{i};
    
    % Separate the RGB channels: each will have dimensions m x n x k
    R = double(currentStack(:,:,1,:));  % Lpc channel (m x n x k)
    G = double(currentStack(:,:,2,:));  % Thymc channel (m x n x k)
    B = double(currentStack(:,:,3,:));  % Blue channel (m x n x k)
    
    % Perform the subtractions for each slice:
    % 1. Subtract the R stack from the G stack (G - R)
    GR_diff = G - R / control_correction_factors(1); 
    
    % 2. Subtract the G stack from the R stack (R - G)
    RG_diff = R - G * control_correction_factors(1);
    
    % 3. Subtract both the R and G stacks from the B stack (B - (R + G))
    BG_diff = B - (R / control_correction_factors(2) + G / control_correction_factors(3));
    
    % Recombine the modified channels into an RGB stack
    % The resulting stack has the same number of slices (k) as the original stack
    cleanedStack = cat(3, RG_diff, GR_diff, BG_diff);  % Combine along the 3rd dimension
    
    % Ensure the cleaned stack is the correct 4D shape (m x n x 3 x k)
    % We need to keep the 4th dimension to preserve the slice count
    imageStacks_clean{i} = uint8(cleanedStack);
end

clearvars -except folderNames_old imageStacks_clean

%% Morphological segmentation of cells and vasculature + Computation

% Initialize a structure to hold the results for each image stack region
results = struct;

% Loop through each region in imageStacks_clean
for p = 1:length(imageStacks_clean)
    % Extract the 3D image stack for the current region
    imageStack = imageStacks_clean{p};  % size: m x n x 3 x k
    
    % Define thresholds for each channel (adjust these thresholds based on your data)
    lowerThresholdLpc = 2;  % Minimum intensity for Lpc population
    upperThresholdLpc = 255;  % Maximum intensity for Lpc population

    lowerThresholdThymc = 3;  % Minimum intensity for Thymc population
    upperThresholdThymc = 255;  % Maximum intensity for Thymc population

    lowerThresholdBlue = 2;  % Minimum intensity for Thymc population
    upperThresholdBlue = 255;  % Maximum intensity for Thymc population

    % Initialize 3D masks for Lpc and Thymc populations (same size as the 3D volume)
    LpcMask3D = false(size(imageStack, 1), size(imageStack, 2), size(imageStack, 4));
    ThymcMask3D = false(size(imageStack, 1), size(imageStack, 2), size(imageStack, 4));
    blueMask3D = false(size(imageStack, 1), size(imageStack, 2), size(imageStack, 4));

    % Loop over all slices in the 3D stack and create 3D masks
    for i = 1:size(imageStack, 4)
        img = imageStack(:,:,:,i);  % Extract the i-th slice
        
        % Segment the populations based on the color channels
        LpcMask = (img(:,:,1) >= lowerThresholdLpc) & (img(:,:,1) <= upperThresholdLpc);
        ThymcMask = (img(:,:,2) >= lowerThresholdThymc) & (img(:,:,2) <= upperThresholdThymc);
        blueMask = (img(:,:,3) >= lowerThresholdBlue) & (img(:,:,3) <= upperThresholdBlue);
        
        % Store the 2D masks into the 3D mask arrays (one for Lpc, one for Thymc)
        LpcMask3D(:,:,i) = LpcMask;
        ThymcMask3D(:,:,i) = ThymcMask;
        blueMask3D(:,:,i) = blueMask;
    end

    % Clean the 3D masks using morphological operations (optional)
    LpcMask3D = imopen(LpcMask3D, strel('disk', 2)); % Clean Lpc mask
    ThymcMask3D = imopen(ThymcMask3D, strel('disk', 4)); % Clean Thymc mask
    blueMask3D = imopen(blueMask3D, strel('disk', 2)); % Clean Thymc mask
    blueMask3D = imclose(blueMask3D, strel('disk', 5)); % Clean Thymc mask
    blueMask3D = imfill(blueMask3D, 'holes');  % Fill holes in blue mask


    % Label the connected components across the 3D volume for both populations
    LpcLabels = bwlabeln(LpcMask3D);
    ThymcLabels = bwlabeln(ThymcMask3D);
    blueLabels = bwlabeln(blueMask3D);
    blueLabels(blueLabels>0) = 1;

    % Count the number of cells in the Lpc and Thymc populations
    numLpcCells = max(LpcLabels(:));  % Number of Lpc cells
    numThymcCells = max(ThymcLabels(:));  % Number of Thymc cells
    

    % Extract cell properties (e.g., centroids) from the 3D masks
    LpcStats = regionprops3(LpcLabels, 'Centroid');
    ThymcStats = regionprops3(ThymcLabels, 'Centroid');

    % Extract centroids for both populations
    LpcCentroids = LpcStats.Centroid; % 3D centroids of Lpc cells
    ThymcCentroids = ThymcStats.Centroid; % 3D centroids of Thymc cells

    results(p).VascularMask = blueLabels;

    % Define scaling factors for each dimension (in microns per pixel)
    scaleXY = 290 / 380;  % Scaling factor for x and y (in microns per pixel)
    scaleZ = 3;           % Scaling factor for z (in microns per pixel)

    % Get the size of the original image stack
    [height, width, ~, numSlices] = size(imageStack);

    % Calculate the physical dimensions of the image in microns
    X_length = width * scaleXY;  % X dimension in microns
    Y_length = height * scaleXY; % Y dimension in microns
    Z_length = numSlices * scaleZ;  % Z dimension in microns
    
    % Calculate the total volume (in cubic microns)
    imagedVolume = X_length * Y_length * Z_length;

    results(p).scaleXY = scaleXY;
    results(p).scaleZ = scaleZ;

    results(p).XdimNative = width;
    results(p).YdimNative = height;
    results(p).ZdimNative = numSlices;

    results(p).XdimInUM = X_length;
    results(p).YdimInUM = Y_length;
    results(p).ZdimInUM = Z_length;

    % Store results in the structure
    results(p).totalVolume = imagedVolume;
    results(p).numLpcCells = numLpcCells;
    results(p).numThymcCells = numThymcCells;

    results(p).LpcCellDensity = results(p).numLpcCells / imagedVolume;
    results(p).ThymcCellDensity = results(p).numThymcCells / imagedVolume;



    % Check if there are no Lpc cells
    if numLpcCells == 0
        results(p).LpcCentroids = [];
    else
        results(p).LpcCentroids = LpcCentroids;  % Scaling Lpc centroids
    end

    % Check if there are no Thymc cells
    if numThymcCells == 0
        results(p).ThymcCentroids = [];
    else
        results(p).ThymcCentroids = ThymcCentroids;  % Scaling Thymc centroids
    end

    % Check if there are no Lpc cells
    if numLpcCells == 0
        results(p).scaledLpcCentroids = [];
    else
        results(p).scaledLpcCentroids = LpcCentroids .* [scaleXY, scaleXY, scaleZ];  % Scaling Lpc centroids
    end

    % Check if there are no Thymc cells
    if numThymcCells == 0
        results(p).scaledThymcCentroids = [];
    else
        results(p).scaledThymcCentroids = ThymcCentroids .* [scaleXY, scaleXY, scaleZ];  % Scaling Thymc centroids
    end
    
end


% Initialize new cell array for folder names
folderNames = cell(size(folderNames_old));

% Iterate through the original cell array
for i = 1:length(folderNames_old)
    % Extract the prefix (e.g., "Tx1", "Tx2", etc.)
    prefix = regexp(folderNames_old{i}, '^[A-Za-z]+\d+', 'match', 'once');
    folderNames{i} = prefix(1:2);
end

% Assuming results is your struct with multiple rows
numRows = length(results);  % Get the number of rows (entries)

% Get the field names from the struct (same for all rows)
fieldNames = fieldnames(results);

% Initialize a cell array to store the data for all rows, plus 1 extra row for the field names
results_cell = cell(numRows + 1, length(fieldNames) + 2);  % Adding space for folderNames and folderNames_old

% Add the "folderNames" and "folderNames_old" as the first two columns in the header row
results_cell(1, 1) = {'folderNames'};  % Header for folderNames
results_cell(1, 2) = {'folderNames_old'};  % Header for folderNames_old

% Add the field names as the remaining columns starting from the 3rd column
results_cell(1, 3:end) = fieldNames';  % Add field names starting from the 3rd column

% Loop through each row of the struct and convert it to a cell array
for row = 1:numRows
    % Add folder names as the first two columns in the results_cell
    results_cell{row + 1, 1} = folderNames{row};
    results_cell{row + 1, 2} = folderNames_old{row};

    % Loop through the field names and store data in corresponding cells
    for fieldIdx = 1:length(fieldNames)
        % Store the data for each field in the corresponding cell
        results_cell{row + 1, fieldIdx + 2} = results(row).(fieldNames{fieldIdx});
    end
end

% Loop through the first column of the results_cell array
for i = 1:size(results_cell, 1)
    if strcmp(results_cell{i, 1}, 'Ct')
        results_cell{i, 1} = 'Control';
    elseif strcmp(results_cell{i, 1}, 'Tk')
        results_cell{i, 1} = 'Sham';
    elseif strcmp(results_cell{i, 1}, 'Tx')
        results_cell{i, 1} = 'Thymc';
    end
end

clearvars -except folderNames_old imageStacks_clean folderNames results_cell

%% Monte carlo nearest neighbor distance and volume overlap calculations


% Define the ellipsoid volume function
ellipsoidVolume = @(a, b, c) (4/3) * pi * a * b * c;

% Number of random samples for volume approximation and number of iterations to use in Monte Carlo integration
numSamples = 1e6; 
numIterations = 1e4;

% Initialize distance_array to store results
distance_array = cell(size(results_cell, 1), 1);  % Cell array to store results for each region

% Identify the column index for 'ThymcCellDensity' (make sure it exists in the header)
ThymcCellDensityColumn = find(strcmp(results_cell(1, :), 'numThymcCells'));

% Step 2: Iterate over each region in results_cell and perform analysis only for Thymc

% Starting from 2 assuming the first row contains headers
for p = 2:size(results_cell, 1)  
    % Check if the region belongs to "Thymc" group
    if strcmp(results_cell{p, 1}, 'Thymc')
        
        % Extract the ThymcCellDensity for the current region
        ThymcCellDensity = cell2mat(results_cell(p, ThymcCellDensityColumn));
        
        % Check if the ThymcCellDensity is above the calculated average
        if ThymcCellDensity > 1 %averageDensity
            
            % Extract Lpc and Thymc Centroids
            LpcCentroids = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'scaledLpcCentroids'))));
            ThymcCentroids = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'scaledThymcCentroids'))));

            % Compute the centroid of Lpc centroids
            LpcCentroidOfCentroidsMean = mean(LpcCentroids, 1);  % Average across rows (X, Y, Z coordinates)
            LpcCentroidOfCentroidsStd = std(LpcCentroids, 1);

            % Compute the centroid of Thymc centroids
            ThymcCentroidOfCentroidsMean = mean(ThymcCentroids, 1);  % Average across rows (X, Y, Z coordinates)
            ThymcCentroidOfCentroidsStd = std(ThymcCentroids, 1);

            % ellipsoid and intersection
            % Extract centroid and radii for the first ellipsoid (Lpc)
            centroid1 = LpcCentroidOfCentroidsMean;  % [x0, y0, z0]
            radii1 = LpcCentroidOfCentroidsStd;     % [a, b, c]
            
            % Extract centroid and radii for the second ellipsoid (Thymc)
            centroid2 = ThymcCentroidOfCentroidsMean;  % [x0, y0, z0]
            radii2 = ThymcCentroidOfCentroidsStd;     % [a, b, c]
        
            % Calculate the volume of each ellipsoid
            volume1 = ellipsoidVolume(radii1(1), radii1(2), radii1(3));
            volume2 = ellipsoidVolume(radii2(1), radii2(2), radii2(3));
        
            % Define the bounds of the bounding box (extending to include both ellipsoids)
            xRange = [min(centroid1(1) - radii1(1), centroid2(1) - radii2(1)), max(centroid1(1) + radii1(1), centroid2(1) + radii2(1))];
            yRange = [min(centroid1(2) - radii1(2), centroid2(2) - radii2(2)), max(centroid1(2) + radii1(2), centroid2(2) + radii2(2))];
            zRange = [min(centroid1(3) - radii1(3), centroid2(3) - radii2(3)), max(centroid1(3) + radii1(3), centroid2(3) + radii2(3))];
        
            % Generate random points within the bounding box
            xRand = xRange(1) + (xRange(2) - xRange(1)) * rand(numSamples, 1);
            yRand = yRange(1) + (yRange(2) - yRange(1)) * rand(numSamples, 1);
            zRand = zRange(1) + (zRange(2) - zRange(1)) * rand(numSamples, 1);
        
            % Define the ellipsoid equations for inclusion test
            isInEllipsoid1 = ((xRand - centroid1(1)).^2 / radii1(1)^2 + (yRand - centroid1(2)).^2 / radii1(2)^2 + (zRand - centroid1(3)).^2 / radii1(3)^2) <= 1;
        
            isInEllipsoid2 = ((xRand - centroid2(1)).^2 / radii2(1)^2 + (yRand - centroid2(2)).^2 / radii2(2)^2 + (zRand - centroid2(3)).^2 / radii2(3)^2) <= 1;
        
            % Calculate the intersection
            intersectionPoints = isInEllipsoid1 & isInEllipsoid2;
            boundedPoints = isInEllipsoid1 | isInEllipsoid2;
        
            % Approximate the fraction of points that are in both ellipsoids
            intersectionFraction = sum(intersectionPoints) / numSamples;
            totalVolumeFraction = sum(boundedPoints) / numSamples;
        
            % Approximate the intersection volume
            intersectionVolume = intersectionFraction * (xRange(2) - xRange(1)) * (yRange(2) - yRange(1)) * (zRange(2) - zRange(1));
            totalboundedVolume = totalVolumeFraction * (xRange(2) - xRange(1)) * (yRange(2) - yRange(1)) * (zRange(2) - zRange(1));
            % intersectionVolume = sum(intersectionPoints) / sum(boundedPoints);
        
            % % Fraction of the volumes that overlap
            fractionOverlap = intersectionVolume / totalboundedVolume;
            
            % Get the bounds (Xdim, Ydim, Zdim)
            XdimInUM = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'XdimInUM'))));
            YdimInUM = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'YdimInUM'))));
            ZdimInUM = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'ZdimInUM'))));

            % Extract the Vascular Mask
            vascularMask = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'VascularMask'))));

            % Match the axial and lateral sampling resolution
            upsampledVascularMask = imresize3(vascularMask, [size(vascularMask, 1), size(vascularMask, 2), 4 * size(vascularMask, 3)]);
            
            % Find the extravascular space
            extravascularIndices = find(upsampledVascularMask == 0);

            % Initialize an array to store the average distances for each iterationred
            avgThymcToLpcDistances = zeros(numIterations, 1);
            avgLpcToThymcDistances = zeros(numIterations, 1);
            fractionOverlap1_all = zeros(numIterations, 1);
            fractionOverlap2_all = zeros(numIterations, 1);
            fractionOverlapMC_all = zeros(numIterations, 1);

            % Check if centroids exist for each population
            if ~isempty(LpcCentroids) && ~isempty(ThymcCentroids)

                % Calculate Thymc to nearest Lpc distance
                ThymcToLpcDistances = zeros(size(ThymcCentroids, 1), 1);
                for j = 1:size(ThymcCentroids, 1)
                    dist = sqrt(sum((LpcCentroids - ThymcCentroids(j, :)).^2, 2));  % Euclidean distance to each Lpc centroid
                    ThymcToLpcDistances(j) = min(dist);  % Minimum distance to a Lpc centroid
                end
                avgThymcToLpcDistance = mean(ThymcToLpcDistances);

                % Calculate Lpc to nearest Thymc distance
                LpcToThymcDistances = zeros(size(LpcCentroids, 1), 1);
                for j = 1:size(LpcCentroids, 1)
                    dist = sqrt(sum((ThymcCentroids - LpcCentroids(j, :)).^2, 2));  % Euclidean distance to each Lpc centroid
                    LpcToThymcDistances(j) = min(dist);  % Minimum distance to a Lpc centroid
                end
                avgLpcToThymcDistance = mean(LpcToThymcDistances);

                % Calculate Thymc to Thymc distance
                ThymcToThymcDistances = zeros(size(ThymcCentroids, 1), 1);
                for i = 1:size(ThymcCentroids, 1)
                    dist = sqrt(sum((ThymcCentroids - ThymcCentroids(i, :)).^2, 2));  % Euclidean distance to all other centroids
                    dist(i) = inf;  % Set the distance to itself to infinity
                    ThymcToThymcDistances(i) = min(dist);  % Nearest neighbor
                end
                avgThymcToThymcDistance = mean(ThymcToThymcDistances);

                % Calculate Lpc to Lpc distance
                LpcToLpcDistances = zeros(size(LpcCentroids, 1), 1);
                for i = 1:size(LpcCentroids, 1)
                    dist = sqrt(sum((LpcCentroids - LpcCentroids(i, :)).^2, 2));  % Euclidean distance to all other centroids
                    dist(i) = inf;  % Set the distance to itself to infinity
                    LpcToLpcDistances(i) = min(dist);  % Nearest neighbor
                end
                avgLpcToLpcDistance = mean(LpcToLpcDistances);

                % Randomize Lpc Centroids and compute Monte Carlo distances
                for i = 1:numIterations
                    % Randomly pick positions within valid regions defined by the upsampled vascular mask
                    randomIndices = extravascularIndices(randi(length(extravascularIndices), size(LpcCentroids, 1), 1));
                    [randomX, randomY, randomZ] = ind2sub(size(upsampledVascularMask), randomIndices);
                    
                    randomizedLpcCentroids = [randomX * (XdimInUM / size(vascularMask, 1)), ...
                                              randomY * (YdimInUM / size(vascularMask, 2)), ...
                                              randomZ * (ZdimInUM / (4 * size(vascularMask, 3)))];
                    
                    % Compute the centroid of randomized Lpc centroids
                    randomLpcCentroid = mean(randomizedLpcCentroids, 1);
                    randomLpcCentroidSTD = std(randomizedLpcCentroids, 1);

                    distances = zeros(size(ThymcCentroids, 1), 1);
                    for j = 1:size(ThymcCentroids, 1)
                        dist = sqrt(sum((randomizedLpcCentroids - ThymcCentroids(j, :)).^2, 2));
                        distances(j) = min(dist);
                    end
                    avgThymcToLpcDistances(i) = mean(distances);

                    distances = zeros(size(randomizedLpcCentroids, 1), 1);
                    for j = 1:size(randomizedLpcCentroids, 1)
                        dist = sqrt(sum((randomizedLpcCentroids(j, :) - ThymcCentroids).^2, 2));
                        distances(j) = min(dist);
                    end
                    avgLpcToThymcDistances(i) = mean(distances);

                    % ellipsoid and intersection
                    % Extract centroid and radii for the first ellipsoid (Lpc)
                    centroid1 = randomLpcCentroid;  % [x0, y0, z0]
                    radii1 = randomLpcCentroidSTD;     % [a, b, c]
                    
                    % Extract centroid and radii for the second ellipsoid (Thymc)
                    centroid2 = ThymcCentroidOfCentroidsMean;  % [x0, y0, z0]
                    radii2 = ThymcCentroidOfCentroidsStd;     % [a, b, c]
                
                    % Calculate the volume of each ellipsoid
                    volume1 = ellipsoidVolume(radii1(1), radii1(2), radii1(3));
                    volume2 = ellipsoidVolume(radii2(1), radii2(2), radii2(3));
                
                    % Define the bounds of the bounding box (extending to include both ellipsoids)
                    xRange = [min(centroid1(1) - radii1(1), centroid2(1) - radii2(1)), max(centroid1(1) + radii1(1), centroid2(1) + radii2(1))];
                    yRange = [min(centroid1(2) - radii1(2), centroid2(2) - radii2(2)), max(centroid1(2) + radii1(2), centroid2(2) + radii2(2))];
                    zRange = [min(centroid1(3) - radii1(3), centroid2(3) - radii2(3)), max(centroid1(3) + radii1(3), centroid2(3) + radii2(3))];
                
                    % Generate random points within the bounding box
                    xRand = xRange(1) + (xRange(2) - xRange(1)) * rand(numSamples, 1);
                    yRand = yRange(1) + (yRange(2) - yRange(1)) * rand(numSamples, 1);
                    zRand = zRange(1) + (zRange(2) - zRange(1)) * rand(numSamples, 1);
                
                    % Define the ellipsoid equations for inclusion test
                    isInEllipsoid1 = ((xRand - centroid1(1)).^2 / radii1(1)^2 + (yRand - centroid1(2)).^2 / radii1(2)^2 + (zRand - centroid1(3)).^2 / radii1(3)^2) <= 1;
                
                    isInEllipsoid2 = ((xRand - centroid2(1)).^2 / radii2(1)^2 + (yRand - centroid2(2)).^2 / radii2(2)^2 + (zRand - centroid2(3)).^2 / radii2(3)^2) <= 1;
                
                    % Calculate the intersection
                    intersectionPoints = isInEllipsoid1 & isInEllipsoid2;
                    boundedPoints = isInEllipsoid1 | isInEllipsoid2;
                
                    % Approximate the fraction of points that are in both ellipsoids
                    intersectionFraction = sum(intersectionPoints) / numSamples;
                    totalVolumeFraction = sum(boundedPoints) / numSamples;
                
                    % Approximate the intersection volume
                    intersectionVolume = intersectionFraction * (xRange(2) - xRange(1)) * (yRange(2) - yRange(1)) * (zRange(2) - zRange(1));
                    totalboundedVolume = totalVolumeFraction * (xRange(2) - xRange(1)) * (yRange(2) - yRange(1)) * (zRange(2) - zRange(1));
                    % intersectionVolume = sum(intersectionPoints) / sum(boundedPoints);
                
                    % % Fraction of the volumes that overlap
                    fractionOverlapMC = intersectionVolume / totalboundedVolume;

                    fractionOverlapMC_all(i) = fractionOverlapMC;

                end
                avgThymcToLpcMonteCarloDist = mean(avgThymcToLpcDistances);
                stdThymcToLpcMonteCarloDist = std(avgThymcToLpcDistances);
                avgLpcToThymcMonteCarloDist = mean(avgLpcToThymcDistances);
                stdLpcToThymcMonteCarloDist = std(avgLpcToThymcDistances);

                avgMonteCarloIntersectMC = mean(fractionOverlapMC_all);

                stdMonteCarloIntersectMC = std(fractionOverlapMC_all);

                % Store results in distance_array
                distance_array{p, 1} = {avgLpcToLpcDistance, avgThymcToThymcDistance, avgThymcToLpcDistance, avgLpcToThymcDistance, ...
                    avgThymcToLpcMonteCarloDist, stdThymcToLpcMonteCarloDist, avgLpcToThymcMonteCarloDist, stdLpcToThymcMonteCarloDist, ...
                    LpcCentroidOfCentroidsMean, LpcCentroidOfCentroidsStd, ThymcCentroidOfCentroidsMean, ThymcCentroidOfCentroidsStd, ...
                    fractionOverlap, avgMonteCarloIntersectMC, stdMonteCarloIntersectMC};
            end
        end
    end
end

% Step 3: Convert to numerical array

% Initialize a matrix to hold the results
numRegions = size(distance_array, 1);  % Number of regions
resultsMatrix = [];  % Empty matrix to store the results

% Loop through the distance_array to extract the values
for p = 2:numRegions  % Start from 2 assuming row 1 contains headers
    if ~isempty(distance_array{p, 1})  % Check if there are results for this region
        % Extract the values from the distance_array (measurements)
        avgLpcToLpcDistance = distance_array{p, 1}{1};
        avgThymcToThymcDistance = distance_array{p, 1}{2};
        avgThymcToLpcDistance = distance_array{p, 1}{3};
        avgLpcToThymcDistance = distance_array{p, 1}{4};
        avgThymcToLpcMonteCarloDist = distance_array{p, 1}{5};
        stdThymcToLpcMonteCarloDist = distance_array{p, 1}{6};
        avgLpcToThymcMonteCarloDist = distance_array{p, 1}{7};
        stdLpcToThymcMonteCarloDist = distance_array{p, 1}{8};
        LpcCentroidOfCentroidsMean = distance_array{p, 1}{9};
        LpcCentroidOfCentroidsStd = distance_array{p, 1}{10};
        ThymcCentroidOfCentroidsMean = distance_array{p, 1}{11};
        ThymcCentroidOfCentroidsStd = distance_array{p, 1}{12};
        fractionOverlap = distance_array{p, 1}{13};
        avgMonteCarloIntersectMC = distance_array{p, 1}{14};
        stdMonteCarloIntersectMC = distance_array{p, 1}{15};

        % Create a row with the p index and the measurements
        resultRow = [p, avgLpcToLpcDistance, avgThymcToThymcDistance, avgThymcToLpcDistance, avgLpcToThymcDistance, ...
                    avgThymcToLpcMonteCarloDist, stdThymcToLpcMonteCarloDist, avgLpcToThymcMonteCarloDist, stdLpcToThymcMonteCarloDist, ...
                    LpcCentroidOfCentroidsMean, LpcCentroidOfCentroidsStd, ThymcCentroidOfCentroidsMean, ThymcCentroidOfCentroidsStd, ...
                    fractionOverlap, avgMonteCarloIntersectMC, stdMonteCarloIntersectMC];
        
        % Append the resultRow to the resultsMatrix
        resultsMatrix = [resultsMatrix; resultRow];
    end
end

% Create a cell array of the nearest neighbor and volume overlap results
headers = {'RegionID', ...
           'LPC_to_Nearest_LPC_Distance', ...
           'ThyMC_to_Nearest_ThyMC_Distance', ...
           'ThyMC_to_Nearest_LPC_Distance', ...
           'LPC_to_Nearest_ThyMC_Distance', ...
           'Average_ThyMC_to_Nearest_RandomDot_Distance', ...
           'Spread_ThyMC_to_Nearest_RandomDot_Distance', ...
           'Average_RandomDot_to_Nearest_ThyMC_Distance', ...
           'Spread_RandomDot_to_Nearest_ThyMC_Distance', ...
           'LPC_DistributionCentroid_X', ...
           'LPC_DistributionCentroid_Y', ...
           'LPC_DistributionCentroid_Z', ...
           'LPC_DistributionSpread_X', ...
           'LPC_DistributionSpread_Y', ...
           'LPC_DistributionSpread_Z', ...
           'ThyMC_DistributionCentroid_X', ...
           'ThyMC_DistributionCentroid_Y', ...
           'ThyMC_DistributionCentroid_Z', ...
           'ThyMC_DistributionSpread_X', ...
           'ThyMC_DistributionSpread_Y', ...
           'ThyMC_DistributionSpread_Z', ...
           'LPC_to_ThyMC_VolumeOverlap', ...
           'Average_RandomDot_to_ThyMC_VolumeOverlap', ...
           'Spread_RandomDot_to_ThyMC_VolumeOverlap'};
       
results_Dist_Overlap_MonteCarlo = cell(1, length(headers));
results_Dist_Overlap_MonteCarlo(1, :) = headers;

results_Dist_Overlap = vertcat(results_Dist_Overlap_MonteCarlo,num2cell(resultsMatrix));

clearvars -except imageStacks_clean results_cell results_Dist_Overlap

%% Plotting the distribution ellipsoids and segmented cell centroids

% Select a region from results_cell
p=21;

% Extract centroid and radii for the ellipsoids
centroid1Plot_X = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'LPC_DistributionCentroid_X'))));
centroid1Plot_Y = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'LPC_DistributionCentroid_Y'))));
centroid1Plot_Z = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'LPC_DistributionCentroid_Z'))));

centroid2Plot_X = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'ThyMC_DistributionCentroid_X'))));
centroid2Plot_Y = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'ThyMC_DistributionCentroid_Y'))));
centroid2Plot_Z = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'ThyMC_DistributionCentroid_Z'))));

radii1Plot_X = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'LPC_DistributionSpread_X'))));
radii1Plot_Y = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'LPC_DistributionSpread_Y'))));
radii1Plot_Z = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'LPC_DistributionSpread_Z'))));

radii2Plot_X = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'ThyMC_DistributionSpread_X'))));
radii2Plot_Y = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'ThyMC_DistributionSpread_Y'))));
radii2Plot_Z = cell2mat(results_Dist_Overlap(p-19, find(strcmp(results_Dist_Overlap(1, :), 'ThyMC_DistributionSpread_Z'))));


LpcCentroidsPlot = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'scaledLpcCentroids'))));
ThymcCentroidsPlot = cell2mat(results_cell(p, find(strcmp(results_cell(1, :), 'scaledThymcCentroids'))));

% Create a grid of points for the first ellipsoid
[x1, y1, z1] = ellipsoid(centroid1Plot_X, centroid1Plot_Y, centroid1Plot_Z, radii1Plot_X, radii1Plot_Y, radii1Plot_Z, 50);

% Create a grid of points for the second ellipsoid
[x2, y2, z2] = ellipsoid(centroid2Plot_X, centroid2Plot_Y, centroid2Plot_Z, radii2Plot_X, radii2Plot_Y, radii2Plot_Z, 50);

% Create the figure
figure;

% Plot the first ellipsoid (Lpc) - Set face color and make sure it doesn't interpolate
h1 = surf(x1, y1, z1);
set(h1, 'FaceColor', 'r', 'EdgeColor', 'none');   % Set color to Lpc (no edge color)

hold on;

% Plot the second ellipsoid (Thymc) - Set face color and make sure it doesn't interpolate
h2 = surf(x2, y2, z2);
set(h2, 'FaceColor', 'g', 'EdgeColor', 'none');   % Set color to Thymc (no edge color)

% Customize the plot
% shading interp;        % Smooth shading for better appearance
% axis equal;            % Equal scaling on all axes
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Two Ellipsoids: Lpc and Thymc');

% Optional: Adjust transparency for better visualization
set(h1, 'FaceAlpha', 0.4);  % Set transparency for the first ellipsoid (Lpc)
set(h2, 'FaceAlpha', 0.4);  % Set transparency for the second ellipsoid (Thymc)

% Plot Thymc Centroids in Thymc color (original ones)
scatter3(ThymcCentroidsPlot(:,1), ThymcCentroidsPlot(:,2), ThymcCentroidsPlot(:,3), ...
         50, 'g', 'filled', 'MarkerFaceAlpha', 0.5);  % 50 is marker size, 'g' is color for Thymc
hold on;

% Plot the original Lpc Centroids in Lpc color (before randomization)
scatter3(LpcCentroidsPlot(:,1), LpcCentroidsPlot(:,2), LpcCentroidsPlot(:,3), ...
         50, 'r', 'filled', 'MarkerFaceAlpha', 0.5);  % 50 is marker size, 'r' is color for Lpc

% Add labels and title
xlabel('X Position (µm)', 'Color', 'k');  % White color for labels
ylabel('Y Position (µm)', 'Color', 'k');
zlabel('Depth (µm)', 'Color', 'k');
set(gca, 'ZDir', 'reverse');
set(gcf, 'Color', 'k');       % Set figure background to black
set(gca, 'Color', 'k');       % Set axes background to black
set(gca, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');  % Axes lines and labels in white
% Set the view to 3D
view(3);
grid on;  % Optional: Add grid for better visibility


clearvars -except imageStacks_clean results_cell results_Dist_Overlap
