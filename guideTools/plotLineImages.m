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
function plotLineImages(segmentsAcum,n,calib,modelClass,u_0,v_0)

theta = -pi:pi/360:pi;

for k = 1:numel(segmentsAcum);
    xHat = segmentsAcum{k};
    plot(xHat(1,:)+u_0,xHat(2,:)+v_0,'r.');
    xHatParam=getParametricCurve(n(:,k),theta,calib,modelClass);
    
    dd = sqrt((xHatParam(1,2:end)-xHatParam(1,1:end-1)).^2+(xHatParam(2,2:end)-xHatParam(2,1:end-1)).^2);
    iii = find(dd>1e5);    
    if numel(iii)==0        
        plot(xHatParam(1,:)+u_0,xHatParam(2,:)+v_0,'g-','LineWidth',1.5);          
    else
        if(iii(1)>1) && (iii(end)<numel(theta))        
            plot(xHatParam(1,1:(iii(1)-1))+u_0,xHatParam(2,1:(iii(1)-1))+v_0,'g-','LineWidth',1.5);          
            plot(xHatParam(1,(iii(end)+1):end)+u_0,xHatParam(2,(iii(end)+1):end)+v_0,'g-','LineWidth',1.5);          
        
%             thetaIni = theta(iii(1)); thetaEnd = theta(iii(end));
%             xHatParamB=getParametricCurve(n(:,k),thetaIni:pi/36000:thetaEnd,calib,modelClass);
%             plot(xHatParamB(1,:)+u_0,xHatParamB(2,:)+v_0,'g.','MarkerSize',5);                    
        else
            thetaIni = 1; thetaEnd = numel(theta);
            xHatParamB=getParametricCurve(n(:,k),thetaIni:pi/36000:thetaEnd,calib,modelClass);
            plot(xHatParamB(1,:)+u_0,xHatParamB(2,:)+v_0,'g.','MarkerSize',5);                                
        end
    end
end
