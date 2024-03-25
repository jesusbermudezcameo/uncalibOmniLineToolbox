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

function [segmentsAcum,nAcum,calib] = uncalibratedExtractionFromBoundaries3PBasedBM(boundaries,subBounds,sizeSubBounds,kEnd,nAttempts,modelClass,lineExtractionParams,rows,u_0,v_0,plotFlag)


shortPixelThreshold = lineExtractionParams.normalizedShortPixelThreshold*rows;

rThreshold = lineExtractionParams.normalizedrThreshold*rows; 


T_center = [[1 0 -u_0];[0 1 -v_0];[0 0 1]];

xAcceptedAcum = cell(1,kEnd);
acceptedIndexMaxAcum = cell(1,kEnd);

r_vlArray = zeros(1,kEnd);
nAcum = zeros(3,kEnd);
fHyper = zeros(1,kEnd);
for j = 1:kEnd
    
    xHat = subBounds(j).xHat;
    gradxF = subBounds(j).gradxF;
    gradyF = subBounds(j).gradyF;
    
    if (plotFlag==1)
        plot(xHat(1,:)+u_0,xHat(2,:)+v_0,'b.','MarkerSize',1);    
    end
        
    if modelClass == 2  
        [acceptedIndex, sizeAccepted, r_vlAcum, fHAcum] = uncalibratedLineExtraction3PointsHyper(xHat,xHat,gradxF,gradyF,nAttempts,shortPixelThreshold,modelClass);        
        fHAcum = fHAcum(~isnan(fHAcum));
        fH = median(fHAcum);
        fHyper(j) = fH;   
        
    else    
        [acceptedIndex, sizeAccepted, r_vlAcum] = uncalibratedLineExtraction3Points(xHat,xHat,nAttempts,shortPixelThreshold,modelClass,rows);
    end
          
    iMax = find(sizeAccepted==max(sizeAccepted));    
    acceptedIndexMax = acceptedIndex{iMax(1)};    
    xAccepted = xHat(:,acceptedIndexMax);             
    
    xAcceptedAcum{j} = xAccepted;
    acceptedIndexMaxAcum{j} = acceptedIndexMax;

    
    r_vlEst = r_vlAcum(iMax(1));
           
    if (size(xAccepted,2)==0)
        errorFlag = 2;
        disp('Unknown Error: uncalibratedExtractionFromBoundaries3PBased.m/72');
    end        
    
    [deltaTheta, phiVar, xSel, theta1, theta2, theta, Rinv, thetaHalf]=getDeltaTheta(xAccepted(1,:),xAccepted(2,:));     

    if (plotFlag==1)
        % Representation
        plot(xAccepted(1,:)+u_0,xAccepted(2,:)+v_0,'r.');    
        plot(xSel(1,:)+u_0,xSel(2,:)+v_0,'b*');          
    end

    if modelClass == 2
        [n, rVL, alpha] = fitCurve(xSel,modelClass,fH); 
    else
        [n, rVL, alpha] = fitCurve(xSel,modelClass); 
    end
    
    
    %% Errors control
    iBreak = 0;
    if ~isreal(alpha) || ~isreal(rVL) || ~isreal(n)
        if(iBreak>5)            
            disp('Assuming bad xSel');
            break;
        end
        disp('bad xSel: uncalibratedExtractionFromBoundaries3PBased.m/71');
        errorFlag = 1;
                
        iRand = ceil(rand(1,2)*size(xAccepted,2));        
        xSelRand = [ xAccepted(:,iRand) xSel(:,3)];                              
        
        if modelClass == 2
            [n, rVL, alpha, l1, l2, l3, r] = fitCurve(xSelRand,modelClass,fH); 
        else
            [n, rVL, alpha, l1, l2, l3, r] = fitCurve(xSelRand,modelClass); %% Ojo que no ira en hiper: ,fHyper);        
        end        
        iBreak = iBreak+1;
    end
    %% End- Errors Control    
    
    r_vlArray(j)=rVL;
             
    
    if size(n,2)==1        
        nAcum(:,j) = n;   
    else
        nAcum(:,j) = [1i;1i;1i];
    end
    
    
end


[r_vlArray, nAcum, xAcceptedAcum, nNan] = avoidNanInf(r_vlArray,nAcum,xAcceptedAcum);
[rVlCand, n, xAccepted, nComplex] = avoidComplex(r_vlArray,nAcum,xAcceptedAcum);

segmentsAcum = xAccepted;

if modelClass == 2
    calib = getCalibFromRVL(median(rVlCand),modelClass,fH);
else
    calib = getCalibFromRVL(median(rVlCand),modelClass);
end


% figure(15);
% plot(rVlCand);
% figure(3);

% deltaThetaSize = zeros(1,numel(xAccepted));
% for k = 1:numel(xAccepted)
%     [deltaTheta,~,subBounds(k).xSelLocal]=getDeltaTheta(xAccepted{k}(1,:),xAccepted{k}(2,:));
%     deltaThetaSize(k)=deltaTheta*sizeSubBounds(k)/rows;
% end

% [deltaThetaSize,iSortKappa]=sort(deltaThetaSize,'descend');
% subBounds = subBounds(iSortKappa);
% sizeSubBounds = sizeSubBounds(iSortKappa);


return;