%     << Automatic Line-Image Extraction Toolbox for Uncalibrated Central
%        Systems with Revolution Symmetry (release v0.5 alpha) >> 
%     ====================================================================
%     Copyright (C) 2014  Jesus Bermudez-Cameo
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%     Acknowledgement:
%     This code has been implemented by Jesus Bermudez-Cameo and supervised 
%     by G.Lopez-Nicolas and J.J. Guerrero
% 
%     One may refer to the following paper:
% 
%     “Line extraction in uncalibrated central images with revolution symmetry”,
%     J. Bermudez-Cameo, G. Lopez-Nicolas and J. J. Guerrero,  
%     24th British Machine Vision Conference, BMVC, pp 1-11; Bristol, UK, Sept.  2013, 
% 
%     This work has been supported by the University of Zaragoza , the 
%     Spanish project VINEA DPI2012-31781, DGA-FSE(T04 and FEDER funds.
%     Jesus Bermudez-Cameo was supported by the FPU program AP2010-3849.

function [boundaries,subBounds,sizeSubBounds,kEnd,boundary] = imagePreprocessingBM(img,mask,preProcessingParams,u_0,v_0)



imgGray=rgb2gray(img);
[rows, columns]=size(imgGray);

cannyIndex = preProcessingParams.cannyIndex;
acumDeltaThetaRatio = preProcessingParams.acumDeltaThetaRatio;  
sizeThresh = preProcessingParams.sizeThresh;
subBoundThresh = preProcessingParams.subBoundThresh; %0.2
filterOmega = preProcessingParams.filterOmega; %2*pi*150

%% Extract edges with canny
[edges, dx, dy] = canny(imgGray,1,[],[1 1]*cannyIndex);
T_center = [[1 0 -u_0];[0 1 -v_0];[0 0 1]];


edges(floor(u_0),floor(v_0))=0;

%% Get connected boundaries and sort by length
[boundaries, boundary, sizeBoundaries]=getBoundaries(edges,mask);
[sortSizeBoundaries, i]=sort(sizeBoundaries,'descend');
boundaries = boundaries(i);


[NUMdd,DENdd] = buildFilter(filterOmega);

subBounds = [];
kSubBounds = 1;
sizeSubBounds = [];
for k = 1:numel(boundaries)
    [binChanges,gradxF,gradyF]= filterBoundSpeedVector(boundaries{k},dx,dy,NUMdd,DENdd,subBoundThresh,rows);    
    [labelSubBounds,nSubBounds]= bwlabel(~binChanges);
    
    for kk = 1:nSubBounds
        subBounds(kSubBounds).xImg = boundaries{k}(:,labelSubBounds==kk);
        subBounds(kSubBounds).gradxF = gradxF(:,labelSubBounds==kk);
        subBounds(kSubBounds).gradyF = gradyF(:,labelSubBounds==kk); 
        sizeSubBounds(kSubBounds) = size(subBounds(kSubBounds).xImg,2);
        kSubBounds = kSubBounds + 1;
    end
    
end

[sizeSubBounds, iSortSubBounds]=sort(sizeSubBounds,'descend');
subBounds = subBounds(iSortSubBounds);
upperI = sizeFilterSegments(sizeSubBounds,sizeThresh(1),sizeThresh(2),sizeThresh(3));
subBounds = subBounds(upperI);
sizeSubBounds = sizeSubBounds(upperI);


%% END COMPUTING THE GRADIENT


deltaThetaSize = zeros(1,numel(subBounds));
for k = 1:numel(subBounds)
    subBounds(k).xHat = T_center*subBounds(k).xImg;
    [deltaTheta,~,subBounds(k).xSelLocal]=getDeltaTheta(subBounds(k).xHat(1,:),subBounds(k).xHat(2,:));
    deltaThetaSize(k)=deltaTheta*sizeSubBounds(k)/rows;
end

[deltaThetaSize,iSortKappa]=sort(deltaThetaSize,'descend');
subBounds = subBounds(iSortKappa);
sizeSubBounds = sizeSubBounds(iSortKappa);

% Select the ones acumulating the 50% of deltaTheta*size

acumSum = sum(deltaThetaSize);
tempSum = 0;
kEnd = 1;
while (tempSum < acumSum*acumDeltaThetaRatio)
    tempSum = sum(deltaThetaSize(1:kEnd));
    kEnd = kEnd + 1;
end




