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

function [deltaTheta phiVar xSel theta1 theta2 theta Rinv thetaHalf]=getDeltaTheta(xHat,yHat)

% normHat = sqrt(xHat.^2+yHat.^2);
% complexvManDir = angle(prod(xHat./normHat+j*yHat./normHat))/numel(xHat);

meanDir = [mean(xHat),mean(yHat)];  %% Plantealo cambiarlo por prod(xHat./normHat+j*yHat./normHat)
normMeanDir = norm(meanDir);  
vManDir = meanDir/normMeanDir;

phiVar = atan2(vManDir(2),vManDir(1));

% R = [[cos(phiVar) -sin(phiVar)];
%      [sin(phiVar) cos(phiVar)]];
 
Rinv = [[cos(phiVar) sin(phiVar)];        %%Rotamos los puntos a una referencia en que el segmento no pase de +pi a -pi consecutivamente.
        [-sin(phiVar) cos(phiVar)]]; 
    
xSeg = Rinv*[xHat;yHat];
theta = atan2(xSeg(2,:),xSeg(1,:));

theta1 = min(theta);
theta2 = max(theta);
deltaTheta = theta2 - theta1;

thetaFromHalf = abs(theta); 
thetaHalf = theta(thetaFromHalf==min(thetaFromHalf));

p = [xHat;yHat;ones(1,numel(xHat))];

x1 = p(:,theta==theta1);
x2 = p(:,theta==theta2);
x3 = p(:,theta==thetaHalf(1));

xSel = [x1(:,1) x2(:,1) x3(:,1)];

return;



