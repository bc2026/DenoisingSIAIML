function JumpFrames = DICAnalysis(Path,ID,X,Y,Sframe,Eframe)
%This function works for a serial DIC analysis. For the sake of simplicity,
%only autotracking is enabled. This program doesn't process movies

%% INPUT variables
%Path: the absolute path of the stream
%ID: the ID name of the partice
%X: X coordinate of particle from ImageJ
%Y: Y coordinate of particle from ImageJ
%Sframe: starting frame of the stream for tracking
%Eframe: ending frame of the stream for tracking

%% OUTPUT
%the number of frames skipped in the tracking
%% settings for analysis (Setup the correct parameter for each run)
% define the particle
endout=regexp(Path,filesep,'split');
stackPath = strjoin(endout(1:end-1),'C:/Users/klee5/OneDrive - stevens.edu/Desktop/Lab Managing/Chloe (SIAI fellowship)');
stackFile = endout{end};
startColX = X+1;% X IN IMAGEJ location of particle in first frame (from ImageJ)
startRowY = Y+1;% Y IN IMAGEJ
particleID = ID;
%Modified 06/15/21 YW
%Adding additional parameters on starting/ending frames for more
%flexibility
SFrame = Sframe; %frame starting for tracking (use 1 as default)
EFrame = Eframe; %frame ending for tracking (use inf as default)

% choose analysis mode
analyzeMode = "track"; % "click","input","track"

% processing parametersok
range = 20; % # of frames to use in moving average. 20 or 6 work well (even)
threshStack = 20; % the threshold for identify none-bkg pixels with all pixels. 40-20, decrease if CD weak, but not below 20 (decrease jumpThresh same time, up to 5)
threshFrame = 0.6; % the threshold fo identify none-bkg pixels within ROI, excluding pixels identified with threshStack.
particleSizeThresh = 20; % # of pixels we consider a particle. cuts out small ones
jumpThresh = 6; % the longest distance allowed in adjacant frames for NP tracking
ROI = 50; % ROI around the particle for final calculation, should be large enough to include movement but not too large to deminish the point of local bkg
ROISmall = 18; % the region around the particle for local bkg calculation
particleDefaultDiameter = 10; % pixels wide of the NP (the size of center circle considered particle)
minArea = floor(pi*(particleDefaultDiameter/2)^2); % pixels (based on circle that is 10 pixel diameter)

% camera parameters
exposureTime=0.100; % unit: s! update this if you are not finding the time stamp from meta data automatically
pixelSize = .065;% um per pixel if 1x zoom using 100x objective

% set up for saving (modify if needed)
[~,saveID,~] = fileparts(stackFile);
savingPath = sprintf('%sProcess/Particle %s/',stackPath,particleID);
if ~isfolder(savingPath) %create folder to save processed data
    mkdir(savingPath);
end
saveBaseName = sprintf('%s%s_Particle%s_range%d_X%d_Y%d',savingPath, saveID, particleID, range, startColX, startRowY);
saveBaseNameMovie = sprintf('%s%s_range%d',savingPath, saveID, range);

%% read image
% since the file processed is not the original file, the time info is lost
% in cropping, so skip the codes for extracting time stamp 
stackName = fullfile(stackPath,stackFile);
stackInfo = imfinfo(stackName); % Getting information from the image 
if ~isfinite(EFrame)
    EFrame = numel(stackInfo);
end
stackFrames = EFrame-SFrame+1; % Getting the number of frames
stackInfo = stackInfo(SFrame:EFrame);%cropping the images
stackMatrix=uint16(zeros(stackInfo(1).Height, stackInfo(1).Width, stackFrames));
for i = 1:stackFrames
    stackMatrix(:,:,i) = imread(stackName,i+SFrame-1);
end

%% subtract moving average background
%This part will average together "range" number of frames before and after
%the frame that is being processed (i.e. a pixel is averaged with the
%same pixel at different points in time).then subtract the average from the
%current frame. This effectively removes conrast from items that are
%CHANGING SLOWER THAN THE PARTICLE, such as the cell and the gradient that
%occurs in DIC images in general.

% pre-allocate matrices
subtractedMatrixDouble=zeros(size(stackMatrix)); %values will be doubles
subtractedMatrix=uint16(zeros(size(stackMatrix))); %values will be unsigned integers
subtractedMatrixB=uint16(zeros(size(stackMatrix)));
subtractedMatrixD=uint16(zeros(size(stackMatrix)));  
for i = 1:stackFrames
    % calculate moving average
    startFrame = max([i-range/2,1]);
    endFrame = min([startFrame+range, stackFrames]);

    if (startFrame - endFrame)<(range-1) %allows for first and last frames to work, but the background subtraction will be different within "range" # of frames
        startFrame = endFrame-range;
    end
    averageImage = uint16(mean(stackMatrix(:,:,startFrame:endFrame),3)); % mean along the 3rd dimension (depth/time)
    %Calculate subtracted image as a double so negative values do not get deleted; aka preserves bright/dark contrast over middle background
    imageMatrixDouble = double(stackMatrix(:,:,i));
    subtractedMatrixDouble(:,:,i) = (imageMatrixDouble-double(averageImage)); %subtract "range" frames from before and after

    subtractedMatrixB(:,:,i)=wiener2(stackMatrix(:,:,i)-averageImage);
    subtractedMatrixD(:,:,i)=wiener2(averageImage-stackMatrix(:,:,i));
    subtractedMatrix(:,:,i)=subtractedMatrixB(:,:,i)+subtractedMatrixD(:,:,i); %this movie will be bright where there are particles and dark everywhere else, for localization purposes, not intensity.

end
minVal=abs(min(min(min(subtractedMatrixDouble)))); % get the darkest pixel in all frames
subtractedMatrixDoubleSaving=uint16(subtractedMatrixDouble+ minVal+100); % make the lowest pixel intensity as 100-the data has to be 'uint16' for saving as image
      

%% Get particle coords - Auto-tracking
% This section only identify the particle with auto-tracking code
%preallocate matrices
particleRegion=false(size(subtractedMatrix)); % label the possible particles (including vesicles) within the whole image
particleRegionROI=false(ROI*2+1, ROI*2+1, stackFrames); % label the possible particles close to the given coords
particleCentroidROI=zeros(stackFrames, 2); %index 1 is X/columns. index 2 is Y/rows.
previousCoord=[ROI+1,ROI+1]; % used to save the identified coords in last frame
fig = figure;
set(fig,'units','normalized','outerposition',[0 0 1 1]);
JumpFrames = 0; %count number of frames with no tracking
for j = 1:stackFrames
    % determine none-bkg pixels from uint16 matrix (all positive int. of
    % NP) under the assumption that most pixels are bkg
    averageUint=mean(mean(subtractedMatrix(:,:,j))); %average of the full frame subtracted matrix (uint16)
    stdevUint=std(std(double(subtractedMatrix(:,:,j)))); 
    % identify particles based on all pixels
    particleRegion(:,:,j)=bwareaopen(subtractedMatrix(:,:,j)>(averageUint+threshStack.*stdevUint), particleSizeThresh); %create logical mask by filtering by intensity and size    
    particleRegionROI(:,:,j)=particleRegion(startRowY-ROI:startRowY+ROI,startColX-ROI:startColX+ROI,j); % possible particles within the ROI
    L=bwlabel(particleRegionROI(:,:,j)); %label the regions in the logical mask so I can select the correct one later if there are more than one

    % get properties of the possible particles regions: the centroid and
    % bounding box
    particleProps = regionprops(L,'Centroid','BoundingBox');

    if j>1
        previousCoord=particleCentroidROI(j-1, :);%coordinate of previous frame as a reference
    end

    % if more than one region identified, choose the one closest to
    % previous centroid (for confined and regular diffusion). For active
    % transport, this can be a problem
    if numel(particleProps)>=2
        coordList = zeros(numel(particleProps),2);
        distances = zeros(numel(particleProps),1);
        for h = 1:numel(particleProps)
%                 coordList(h,:) = [particleProps(h).WeightedCentroid(1); particleProps(h).WeightedCentroid(2)];
            coordList(h,:) = [particleProps(h).Centroid(1); particleProps(h).Centroid(2)];
            distances(h) = pdist([previousCoord; coordList(h,:)]);
            [~,chosenParticle] = min(distances); % chosenParticle saves the index of the coords of smallest distance, under the assumption cellular features are not identified
        end

        particleProps = particleProps(chosenParticle);%set region propes to only include the chosen particle
        particleRegionROI(:,:,j) = L==chosenParticle;%set particle mask  to only include the chosen particle
    end %numel(particleProps) should be 1 or 0 now after this if statement

    if numel(particleProps)==1 % either only find one or picked closest one from last couple lines
%             particleCentroidROI(j,:)=[particleProps(1).WeightedCentroid(1), particleProps(1).WeightedCentroid(2)];
        particleCentroidROI(j,:)=[particleProps(1).Centroid(1), particleProps(1).Centroid(2)];

        % if the particle jumped too far, set the centroid as last frame
        % location and define a 'default particle mask' centered around the
        % centroid.
        distanceFromPrev = pdist([previousCoord; particleCentroidROI(j,:)]);
        if distanceFromPrev>jumpThresh
            JumpFrames = JumpFrames +1;
            particleCentroidROI(j,:)=previousCoord; %assume particle didn't move from last frame
            [rr, cc] = meshgrid(1:2*ROI+1);
            particleRegionROI(:,:,j)=((rr-particleCentroidROI(j,1)).^2+(cc-particleCentroidROI(j,2)).^2)<(particleDefaultDiameter./2).^2;
        end
    end

    if numel(particleProps)==0
        JumpFrames = JumpFrames +1;
        particleCentroidROI(j,:)=previousCoord;
        [rr, cc] = meshgrid(1:2*ROI+1);
        particleRegionROI(:,:,j)=((rr-particleCentroidROI(j,1)).^2+(cc-particleCentroidROI(j,2)).^2)<(particleDefaultDiameter./2).^2;
    end

    % if the area is too small, make it at least min diameter
    particlePropsFinal=regionprops(particleRegionROI(:,:,j));
    if particlePropsFinal.Area<minArea
        [rr, cc] = meshgrid(1:2*ROI+1);
        particleRegionROI(:,:,j)=((rr-particleCentroidROI(j,1)).^2+(cc-particleCentroidROI(j,2)).^2)<(particleDefaultDiameter./2).^2;
    end
end

centroidGlobalColX = round(particleCentroidROI(:,1))+(startColX-(ROI+1));
centroidGlobalRowY = round(particleCentroidROI(:,2))+(startRowY-(ROI+1));

T = table(centroidGlobalColX,centroidGlobalRowY);
writetable(T,sprintf('%s_Coords.xlsx',saveBaseName));


%% Calculation

% preallocate matrices
rawROI=zeros(ROI*2+1, ROI*2+1, stackFrames); % save the small ROI
nonebkgROI=false(ROI*2+1, ROI*2+1, stackFrames); % label the possible particles close to the given coords
bkgregionROI = false(ROI*2+1, ROI*2+1, stackFrames); % label the bkg pixels considered in calculation
localAverage=zeros(stackFrames,1); % local bkg
localStdev=zeros(stackFrames,1); % std for local bkg
particleSignal=zeros(stackFrames,1); % mean intensity of all pixels INSIDE the particle mask
particleSignalBright=zeros(stackFrames,1); % mean intensity of identified bright pixels
particleSignalDark=zeros(stackFrames,1); % mean intensity of identified dark pixels
nanCountBrights=zeros(stackFrames,1); % record whether bright pixels were identified in each frame
nanCountDarks=zeros(stackFrames,1); % record whether dark pixels were identified in each frame
numberBrightPixels=zeros(stackFrames,1); % record # of bright pixels identified in each frame
numberDarkPixels=zeros(stackFrames,1); % record # of dark pixels identified in each frame

for j = 1:stackFrames
    % save the ROI images for later reference
    rawROI(:,:,j)=stackMatrix(startRowY-ROI:startRowY+ROI,startColX-ROI:startColX+ROI,j); % get the bkg subtracted image around initial coords
    rawROIHolder=rawROI(:,:,j); % the current frame
    %subtractedMatrixROIHolder=subtractedMatrix(startRowY-ROI:startRowY+ROI,startColX-ROI:startColX+ROI,j);
    
    %% Locate the high-contrast vesicles
    % re-identify the "particles" in the current frame based on raw image
    % instead of moving average subtracted image (similar to fixed study)
    averageROI = uint16(mean(mean(rawROIHolder)));
    
    % this equals to the subtracted image used before
    NPidentifyROI = medfilt2(uint16(rawROIHolder)-averageROI)+medfilt2(averageROI-uint16(rawROIHolder));
    averageFrame = mean(mean(NPidentifyROI));
    stdevFrame = std(std(double(NPidentifyROI)));
    
    % identify the high-contrast bkg stuff
    nonebkgROI(:,:,j) = bwareaopen(NPidentifyROI>(averageFrame+threshFrame.*stdevFrame), particleSizeThresh); 
    bkgregionROI(:,:,j) = ~nonebkgROI(:,:,j) & ~particleRegionROI(:,:,j); % make sure particle region is also excluded
    
    %% Calculate based on identified particle centroid
    % determine LOCAl bkg intensities in a smaller area
    particleCentroidROIround=round(particleCentroidROI(j,:)); %round centroid location
    if particleCentroidROIround(2)-ROISmall<1 ||particleCentroidROIround(1)-ROISmall<1 ||particleCentroidROIround(2)+ROISmall> ROI*2+1 ||particleCentroidROIround(1)+ ROISmall > ROI*2+1
        disp('Index Error')
        return
    end
    rawROIholderLocal=rawROI(particleCentroidROIround(2)-ROISmall:particleCentroidROIround(2)+ROISmall,particleCentroidROIround(1)-ROISmall:particleCentroidROIround(1)+ROISmall,j);
    notParticlesLocal=bkgregionROI(particleCentroidROIround(2)-ROISmall:particleCentroidROIround(2)+ROISmall,particleCentroidROIround(1)-ROISmall:particleCentroidROIround(1)+ROISmall,j);
    %particleZoomProps = regionprops(particleRegionROI(particleCentroidROIround(2)-ROISmall:particleCentroidROIround(2)+ROISmall,particleCentroidROIround(1)-ROISmall:particleCentroidROIround(1)+ROISmall,j));
    %particleFinalProps=regionprops(particleRegionROI(:,:,j));
    
    localAverage(j)=mean(rawROIholderLocal(notParticlesLocal));%average of local background excluding particle
    localStdev(j)=std(rawROIholderLocal(notParticlesLocal));
    
    % determine mean signal intensity from inside the particle mask. NOTE:
    % use raw image data
    particleSignal(j)=mean(rawROIHolder(particleRegionROI(:,:,j))); %this is the average intensity of all pixels INSIDE the particle mask
    
    particleSignalList=rawROIHolder(particleRegionROI(:,:,j)); %list of intensity values inside the particle mask
    
    particleSignalListBright=particleSignalList>localAverage(j); %logical matrix of pixel values that are greater than background
    particleSignalListDark=particleSignalList<localAverage(j); %logical matrix of pixel values that are less than background
    numberBrightPixels(j)=sum(particleSignalListBright);
    numberDarkPixels(j)=sum(particleSignalListDark);
    
    particleSignalBright(j)=mean(particleSignalList(particleSignalListBright)); %this is the average intensity of pixels inside the particle mask that are ABOVE the local background
    particleSignalDark(j)=mean(particleSignalList(particleSignalListDark));%this is the average intensity of pixels inside the particle mask that are less than the local background

    %if there are no pixels above or below the average (i.e. FULL dark or bright), then replace the signal value the background?.
    if isnan(particleSignalBright(j))
        particleSignalBright(j)=localAverage(j);
        nanCountBrights(j)=1;
    end
    
    if isnan(particleSignalDark(j))
        particleSignalDark(j)=localAverage(j);
        nanCountDarks(j)=1;
    end
    
end

%contrast calculations
brightContrast=(particleSignalBright-localAverage)./localAverage;
darkContrast=(particleSignalDark-localAverage)./localAverage;
bkgContrast=(localAverage-localAverage)./localAverage;

brightContrastNorm=(brightContrast-min(brightContrast))./(max(brightContrast)-min(brightContrast));%normalized bright contrast
darkContrastNorm=(darkContrast-max(darkContrast))./(min(darkContrast)-max(darkContrast));%normalized dark contrast

contrastDifference=brightContrastNorm-darkContrastNorm;

contrast = (particleSignal-localAverage)./localAverage;

%centroid calculations
particleCentroidROIColXSmooth=smoothdata(particleCentroidROI(:,1));%smooth ROI centroids
particleCentroidROIRowYSmooth=smoothdata(particleCentroidROI(:,2));

centroidGlobalColXSmooth=smoothdata(centroidGlobalColX);%smooth global centroids
centroidGlobalRowYSmooth=smoothdata(centroidGlobalRowY);

positionColX_um=centroidGlobalColX.*pixelSize;%X position in um
positionRowY_um=centroidGlobalRowY.*pixelSize;%Y position in um

positionColXFirst=positionColX_um(1);
positionRowYFirst=positionRowY_um(2);

displacement_um=((positionColX_um-positionColXFirst).^2+(positionRowY_um-positionRowYFirst).^2).^.5; %dispalcement from the frist frame in microns
displacement_umSmooth=smoothdata(displacement_um);


%% save the result

frameNum = (1:stackFrames)';
timeStampExp = (0:exposureTime:(stackFrames-1)*exposureTime)';
particleCentroidROIX = particleCentroidROI(:,1);
particleCentroidROIY = particleCentroidROI(:,2);

% save the values for each frame
TFrame = table(frameNum,timeStampExp,localAverage,localStdev,particleSignal,contrast,...
    particleSignalBright,particleSignalDark,brightContrast,darkContrast,brightContrastNorm,darkContrastNorm,contrastDifference,... % contrast difference
    particleCentroidROIX,particleCentroidROIY,particleCentroidROIColXSmooth,particleCentroidROIRowYSmooth,...
    centroidGlobalColX,centroidGlobalRowY,centroidGlobalColXSmooth,centroidGlobalRowYSmooth,... % position of particle
    displacement_um,displacement_umSmooth); %,particleCentroidROIColXSmooth,particleCentroidROIRowYSmooth,centroidGlobalColXSmooth,centroidGlobalRowYSmooth
writetable(TFrame,sprintf('%s_Result.xlsx',saveBaseName))

% save the related parameters for the analysis
imageFolder = string(stackPath);
imageFile = string(stackFile);
TParticle = table(imageFolder,imageFile,SFrame,EFrame,startColX,startRowY,analyzeMode,ROI,ROISmall,range,threshStack,threshFrame,particleSizeThresh,particleDefaultDiameter,minArea,jumpThresh,exposureTime,pixelSize);
writetable(TParticle,sprintf('%s_Parameters.xlsx',saveBaseName))

roistack = uint16(rawROI);
imwrite(roistack(:,:,1),sprintf('%s_ROI.tif',saveBaseName));
for k = 2:size(subtractedMatrixDoubleSaving,3)
    imwrite(roistack(:,:,k), sprintf('%s_ROI.tif',saveBaseName), 'writemode', 'append')
end

particleidentified = zeros(size(roistack));
for h = 1:size(particleidentified,3)
    roiframe = roistack(:,:,h);
    particleframe = zeros(size(roiframe));
    particleframe(particleRegionROI(:,:,h)) = roiframe(particleRegionROI(:,:,h));
    particleframe(~particleRegionROI(:,:,h)) = localAverage(h)-100;
    particleidentified(:,:,h) = particleframe;
end
particleidentified = uint16(particleidentified);
imwrite(particleidentified(:,:,1),sprintf('%s_Particle %s.tif',saveBaseNameMovie,particleID));
for k = 2:size(particleidentified,3)
    imwrite(particleidentified(:,:,k), sprintf('%s_Particle %s.tif',saveBaseNameMovie,particleID), 'writemode', 'append')
end

% if JumpFrames ~=0
%     %fprintf('Total %d frames were not properly tracked in this run.\n',JumpFrames);
% end
%% Plot figures
% signal figures
figure; % plot the average particle signal comparing to local background
plot(timeStampExp,localAverage);
hold on
plot(timeStampExp,localAverage+localStdev);
plot(timeStampExp,localAverage-localStdev)
plot(timeStampExp,particleSignal)
title(sprintf('Mean Signal, file %s, Particle %s, StartX %d, startY %d', saveID, particleID, startColX, startRowY),'interpreter','none')
legend('Local background','+ Local stdev','- Local stdev', 'Particle signal')
xlabel('Time (sec)')
ylabel('Intensity')
% saveas(gcf,sprintf('%s_SignalMean.fig',saveBaseName));
% saveas(gcf,sprintf('%s_SignalMean.png',saveBaseName));
saveas(gcf,sprintf('%s_SignalMean.eps',saveBaseName),'epsc');

figure; % plot the bright and dark pixel of particle comparing to the background
plot(timeStampExp,localAverage);
hold on
plot(timeStampExp,localAverage+localStdev)
plot(timeStampExp,localAverage-localStdev)
plot(timeStampExp,particleSignalBright)
plot(timeStampExp,particleSignalDark)
title(sprintf('Bright and Dark Signal, file %s, Particle %s, startX %d, startY %d',saveID, particleID, startColX, startRowY),'interpreter','none')
legend('Local background','+ Local stdev','- Local stdev','Bright signal','Dark signal')
xlabel('Time (sec)')
ylabel('Intensity')
% saveas(gcf,sprintf('%s_SignalBrightDark.fig',saveBaseName));
% saveas(gcf,sprintf('%s_SignalBrightDark.png',saveBaseName));
saveas(gcf,sprintf('%s_SignalBrightDark.eps',saveBaseName),'epsc');

% contrast figures
figure; % plot the bright and dark contrast
plot(timeStampExp, bkgContrast)
hold on
ax = gca; % these two lines are intended to make bright and dark lines to be consistent for all the plots (purple,green)
ax.ColorOrderIndex = 4;
plot(timeStampExp,brightContrast);
plot(timeStampExp,darkContrast);
title(sprintf('Bright and Dark Contrast, file %s, Particle %s, startX %d, startY %d', saveID, particleID,startColX,startRowY),'interpreter','none')
legend('Background','Bright Contrast','Dark Contrast')
xlabel('Time (sec)')
ylabel('Contrast')
% saveas(gcf,sprintf('%s_BD Contrast.fig',saveBaseName));
% saveas(gcf,sprintf('%s_Contrast.png',saveBaseName));
saveas(gcf,sprintf('%s_BD Contrast.eps',saveBaseName),'epsc');


figure; % plot the normalized bright and dark contrast
hold on
ax = gca;
ax.ColorOrderIndex = 4;
plot(timeStampExp,brightContrastNorm);
plot(timeStampExp,darkContrastNorm);
title(sprintf('Normalized Bright and Dark Contrast, file %s, Particle %s, startX %d, startY %d',saveID,particleID,startColX, startRowY),'interpreter','none')
legend('Normalied Bright Contrast','Normalized Dark Contrast')
xlabel('Time (sec)')
ylabel('Normalized Contrast')
% saveas(gcf,sprintf('%s_ContrastNorm.fig',saveBaseName));
% saveas(gcf,sprintf('%s_ContrastNorm.png',saveBaseName));
saveas(gcf,sprintf('%s_ContrastNorm.eps',saveBaseName),'epsc');


figure; % plot the contrast difference
plot(timeStampExp,contrastDifference);
hold on
plot(timeStampExp,zeros(1,stackFrames)) % line y=0
title(sprintf('Contrast Difference, file %s, Paticle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('Time (sec)')
ylabel('Contrast Difference')
ylim([-1,1])
% saveas(gcf,sprintf('%s_ContrastDiff.fig',saveBaseName));
% saveas(gcf,sprintf('%s_ContrastDiff.png',saveBaseName));
saveas(gcf,sprintf('%s_ContrastDiff.eps',saveBaseName),'epsc');

figure; % plot the contrast difference
plot(timeStampExp,contrast);
hold on
plot(timeStampExp,zeros(1,stackFrames)) % line y=0
title(sprintf('Contrast, file %s, Paticle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('Time (sec)')
ylabel('Contrast')
% saveas(gcf,sprintf('%s_Contrast.fig',saveBaseName));
% saveas(gcf,sprintf('%s_ContrastDiff.png',saveBaseName));
saveas(gcf,sprintf('%s_Contrast.eps',saveBaseName),'epsc');

% displacement figures
figure; % plot the movement in X axis
plot(timeStampExp,centroidGlobalColX)
hold on
plot(timeStampExp,centroidGlobalColXSmooth)
title(sprintf('X position, file %s, Particle %s, startX %d, startY %d', saveID, particleID, startColX, startRowY),'interpreter','none')
xlabel('Time (sec)')
ylabel('X position (pixels)')
legend('x','x-smooth')
% saveas(gcf,sprintf('%s_positionX.fig',saveBaseName));
% saveas(gcf,sprintf('%s_positionX.png',saveBaseName));
saveas(gcf,sprintf('%s_positionX.eps',saveBaseName),'epsc');


figure; % plot the movement in Y axis
plot(timeStampExp,centroidGlobalRowY)
hold on
plot(timeStampExp,centroidGlobalRowYSmooth)
title(sprintf('Y position, file %s, Particle %s, startX %d, startY %d', saveID, particleID, startColX, startRowY),'interpreter','none')
xlabel('Time (sec)')
ylabel('Y position (pixels)')
legend('y','y-smooth')
% saveas(gcf,sprintf('%s_positionY.fig',saveBaseName));
% saveas(gcf,sprintf('%s_positionY.png',saveBaseName));
saveas(gcf,sprintf('%s_positionY.eps',saveBaseName),'epsc');


figure; % plot the displacement of the particle
plot(timeStampExp,displacement_um)
hold on
plot(timeStampExp,displacement_umSmooth)
title(sprintf('Displacement (um), file %s, Particle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('Time (sec)')
ylabel('Displacement (um)')
legend('Displacement','Displacement-smooth')
% saveas(gcf,sprintf('%s_Displacement.fig',saveBaseName));
% saveas(gcf,sprintf('%s_Displacement.png',saveBaseName));
saveas(gcf,sprintf('%s_Displacement.eps',saveBaseName),'epsc');

% overlay figures
figure; % overlay the trajectory with the whole image
imagesc(stackMatrix(:,:,1))
colormap('gray')
axis equal tight
hold on
scatter(centroidGlobalColX(1),centroidGlobalRowY(1)) % the position in first frame
plot(centroidGlobalColX,centroidGlobalRowY) % trajectory
legend('Start position','Path')
title(sprintf('Overlay, file %s, Particle %s, startX %d, startY %d',saveID, particleID, startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_Overlay.fig',saveBaseName));
% saveas(gcf,sprintf('%s_Overlay.png',saveBaseName));
saveas(gcf,sprintf('%s_Overlay.eps',saveBaseName),'epsc');


figure; % overlay smoothed trajectory on whole image
imagesc(stackMatrix(:,:,1));
colormap('gray')
axis equal tight
hold on
scatter(centroidGlobalColX(1),centroidGlobalRowY(1))
plot(centroidGlobalColXSmooth,centroidGlobalRowYSmooth)
legend('Start position','Path')
title(sprintf('Overlay Smooth, file %s, Particle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_OverlaySmooth.fig',saveBaseName));
% saveas(gcf,sprintf('%s_OverlaySmooth.png',saveBaseName));
saveas(gcf,sprintf('%s_OverlaySmooth.eps',saveBaseName),'epsc');


figure; % overlay trajectory in small ROI
imagesc(roistack(:,:,1))
colormap('gray')
axis equal tight
hold on
scatter(particleCentroidROI(1,1),particleCentroidROI(1,2));
plot(particleCentroidROI(:,1),particleCentroidROI(:,2))
legend('Start position','Path')
title(sprintf('Overlay zoom, file %s, Particle %s, startX %d, startY %d', saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_OverlayZoom.fig',saveBaseName));
% saveas(gcf,sprintf('%s_OverlayZoom.png',saveBaseName));
saveas(gcf,sprintf('%s_OverlayZoom.eps',saveBaseName),'epsc');


figure; % overlay smooth trajectory in ROI
imagesc(roistack(:,:,1))
colormap('gray');
axis equal tight
hold on
scatter(particleCentroidROI(1,1),particleCentroidROI(1,2))
plot(particleCentroidROIColXSmooth,particleCentroidROIRowYSmooth)
legend('Start position','Path')
title(sprintf('Overlay zoom smooth, file %s, Particle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_OverlayZoomSmooth.fig',saveBaseName));
% saveas(gcf,sprintf('%s_OverlayZoomSmooth.png',saveBaseName));
saveas(gcf,sprintf('%s_OverlayZoomSmooth.eps',saveBaseName),'epsc');


% plot paths with color according to time, set axes to aesily overlay on
% the raw image or the ROI
figure; % plot the path only
c = 1:length(centroidGlobalColX);
patch([centroidGlobalColX',NaN], [centroidGlobalRowY' NaN], [c NaN], [c NaN], 'edgecolor','interp');
hold on
scatter(centroidGlobalColX(1),centroidGlobalRowY(1))
colormap('jet')
set(gca,'yDir','reverse')
axis equal tight
xlim([0 size(stackMatrix(:,:,1),2)])
ylim([0 size(stackMatrix(:,:,1),1)])
legend('Path','Start position')
title(sprintf('Full path, file %s, Particle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_PathFull.fig',saveBaseName));
% saveas(gcf,sprintf('%s_PathFull.png',saveBaseName));
saveas(gcf,sprintf('%s_PathFull.eps',saveBaseName),'epsc');


figure; % plot path in the ROI
c = 1:length(particleCentroidROI);
patch([particleCentroidROI(:,1)' NaN], [particleCentroidROI(:,2)' NaN], [c NaN], [c NaN], 'edgecolor','interp');
hold on 
scatter(particleCentroidROI(1,1),particleCentroidROI(1,2))
colormap('jet')
set(gca, 'yDir', 'reverse')
axis equal tight
xlim([0 2*ROI+1])
ylim([0 2*ROI+1])
legend('Path','Start position')
title(sprintf('ROI path, file %s, Particle %s, startX %d, startY %d',saveID, particleID, startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_PathROI.fig',saveBaseName));
% saveas(gcf,sprintf('%s_PathROI.png',saveBaseName));
saveas(gcf,sprintf('%s_PathROI.eps',saveBaseName),'epsc');


figure; % plot smooth path only in whole image 
c = 1:length(centroidGlobalColXSmooth);
patch([centroidGlobalColXSmooth' NaN], [centroidGlobalRowYSmooth' NaN], [c NaN], [c NaN],'edgecolor','interp');
hold on
scatter(centroidGlobalColX(1),centroidGlobalRowY(1))
colormap('jet')
set(gca,'yDir','reverse')
axis equal tight
xlim([0 size(stackMatrix(:,:,1),2)])
ylim([0 size(stackMatrix(:,:,1),1)])
legend('Path','Start position')
title(sprintf('Full path smooth, file %s, Particle %s, startX %d, startY %d',saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_PathFullSmooth.fig',saveBaseName));
% saveas(gcf,sprintf('%s_PathFullSmooth.png',saveBaseName));
saveas(gcf,sprintf('%s_PathFullSmooth.eps',saveBaseName),'epsc');


figure; % plot smooth trajectory in ROI
c = 1:length(particleCentroidROIColXSmooth);
patch([particleCentroidROIColXSmooth' NaN], [particleCentroidROIRowYSmooth' NaN],[c NaN],[c NaN],'edgecolor','interp');
hold on
scatter(particleCentroidROI(1,1),particleCentroidROI(1,2))
colormap('jet')
set(gca,'yDir','reverse')
axis equal tight
xlim([0 2*ROI+1])
ylim([0 2*ROI+1])
legend('Path','Start position')
title(sprintf('ROI path smooth, file %s, Particle %s, startX %d, startY %d', saveID, particleID, startColX, startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_PathROISmooth.fig',saveBaseName));
% saveas(gcf,sprintf('%s_PathROISmooth.png',saveBaseName));
saveas(gcf,sprintf('%s_PathROISmooth.eps',saveBaseName),'epsc');


% save the starting frame
figure; % the whole image
imagesc(stackMatrix(:,:,1))
axis equal tight
title(sprintf('Frame %d, file %s, Particle %s, startX %d, startY %d', SFrame,saveID, particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_Frame1.fig',saveBaseName));
% saveas(gcf,sprintf('%s_Frame1.png',saveBaseName));
saveas(gcf,sprintf('%s_Frame%d.eps',saveBaseName,SFrame),'epsc');


figure; % plot the ROI
imagesc(rawROI(:,:,1))
colormap('gray')
axis equal tight
title(sprintf('Raw ROI, file %s, Particle %s, startX %d, startY %d',saveID,particleID,startColX,startRowY),'interpreter','none')
xlabel('X (pixels)')
ylabel('Y (pixels)')
% saveas(gcf,sprintf('%s_Frame1RawROI.fig',saveBaseName));
% saveas(gcf,sprintf('%s_Frame1RawROI.png',saveBaseName));
saveas(gcf,sprintf('%s_Frame1RawROI.eps',saveBaseName),'epsc');

close all
end