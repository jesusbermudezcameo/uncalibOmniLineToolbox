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

function binInliersDist = ransacExtractionFromBoundaryBL(xHatSel,xHatTest,calib,modelClass,rThreshold,nAttempts)

p=2; 
nPointsSel = size(xHatSel,2);
index=ceil(rand(nAttempts,p)*nPointsSel);

rTest = sqrt(xHatTest(1,:).^2+xHatTest(2,:).^2);
alphaTest = getAlpha(rTest,calib,modelClass);

rSel =sqrt(xHatSel(1,:).^2 + xHatSel(2,:).^2);
alphaSel = getAlpha(rSel,calib,modelClass);
s = getSign(modelClass);

xHat1 = xHatSel(:,index(:,1));
xHat2 = xHatSel(:,index(:,2));
alpha1 = alphaSel(index(:,1));
alpha2 = alphaSel(index(:,2));
n = twoPointsToNormalM(xHat1,xHat2,alpha1,alpha2,s);

%% Comprobar
% for k = 1:size(xHat1,2)
%     [n2(:,k),alphaP] = twoPointsToNormal([xHat1(:,k),xHat2(:,k)],calib,modelClass);  
%     alpha1B(k) = alphaP(1);
%     alpha2B(k) = alphaP(2);
% end


dAlg = n'*[xHatTest(1,:);s*xHatTest(2,:);-alphaTest];

dalpha_r = getDalpha_r(rSel(index),calib,modelClass);
[gradx,grady] = gradGeneral(n(1,:),s*n(2,:),n(3,:),[xHat1(1,:);xHat2(1,:)],[xHat1(2,:);xHat2(2,:)],dalpha_r');
algThr = getAlgThreshold([xHat1(1,:);xHat2(1,:)],[xHat1(2,:);xHat2(2,:)],gradx,grady,n,rThreshold,calib,modelClass,s);


binInliersDist = bsxfun(@lt,abs(dAlg),algThr');

return


% dalpha_rTest = getDalpha_r(rTest,calib,modelClass);
% % gradx = nx-nz*dalpha_r*x
% gradxTest = ...
%     bsxfun(@minus,n(1,:)',bsxfun(@times,n(3,:)',repmat(dalpha_rTest.*xHatTest(1,:),size(n,2),1)));
% % grady = ny-nz*dalpha_r*y
% gradyTest = ...
%     bsxfun(@minus,n(2,:)',bsxfun(@times,n(3,:)',repmat(dalpha_rTest.*xHatTest(2,:),size(n,2),1)));
% 
% normGradTest = sqrt(gradxTest.^2+gradyTest.^2);
% gradxTest = gradxTest./normGradTest;
% gradyTest = gradyTest./normGradTest;











function [gradx,grady] = gradGeneral(nx,ny,nz,x,y,dalpha_r)
gradx = bsxfun(@minus,nx,bsxfun(@times,nz,dalpha_r.*x));
grady = bsxfun(@minus,ny,bsxfun(@times,nz,dalpha_r.*y));

normGrad = sqrt(gradx.^2+grady.^2);
gradx = gradx./normGrad;
grady = grady./normGrad;
return
 

function algThr = getAlgThreshold(x,y,gradx,grady,n,rThreshold,calib,modelClass,s)
xEdge = x + gradx*rThreshold;
yEdge = y + grady*rThreshold;
rEdge = sqrt(xEdge.^2+yEdge.^2);
alpha = getAlpha(rEdge,calib,modelClass);  
algThr = mean(abs(bsxfun(@times,n(1,:),xEdge) + s*bsxfun(@times,n(2,:),yEdge) - bsxfun(@times,n(3,:),alpha)));
return
