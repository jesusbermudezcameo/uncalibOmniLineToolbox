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

function  [acceptedIndex, sizeAccepted, r_vlAcum, fHAcum] = uncalibratedLineExtraction3PointsHyper(xHatSel,xHatTest,gradxFSel,gradyFSel,nAttempts,rThreshold,modelClass)

sign=getSign(modelClass);
p=3; 
nPointsSel = size(xHatSel,2);
index=ceil(rand(nAttempts,p)*nPointsSel);

acceptedIndex = cell(1,size(index,1));
sizeAccepted = zeros(1,size(index,1));
r_vlAcum = zeros(1,size(index,1));
fHAcum = zeros(1,size(index,1));

for k=1:size(index,1)
    xRand = xHatSel(:,index(k,:));       
    dxRand = [gradxFSel(index(k,:));gradyFSel(index(k,:))];

    
    fH = gradientsToFocalHyperNoRVL(dxRand,xRand);      
       
    
    [n rVL alpha l1 l2 l3 r] = fitCurve(xRand,modelClass,fH);    
    
    r_vlAcum(k) = rVL;
    fHAcum(k) = fH;
    
    if isreal(rVL)&&~isnan(rVL)&&~isinf(rVL)
    
        calib = getCalibFromRVL(rVL,modelClass,fH); 

        rTest = sqrt(xHatTest(1,:).^2+xHatTest(2,:).^2);    
        %% Computing Distance
        alphaTest = getAlpha(rTest,calib,modelClass);
        dist = abs([xHatTest(1,:);sign*xHatTest(2,:);-alphaTest]'*n);

        %% Computing adaptative Threshold
        grad = getGradient(n,xRand(1,:),xRand(2,:),calib,modelClass);
        dxSel = normalizeVectorArray(grad);
        algThr = getAlgThreshold(xRand,n,dxSel,rThreshold,calib,modelClass,sign);

        iAccepted = find(dist<algThr);
        sizeAccepted(k) = numel(iAccepted);
        acceptedIndex{k} = iAccepted; 
    else
        sizeAccepted(k) = 0;
        acceptedIndex{k} = [];         
    end              
end