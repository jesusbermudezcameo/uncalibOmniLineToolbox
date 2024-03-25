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

function [n, rVL, alpha]=fitCurveHyper(x,r,l1,l2,l3,f)

n = [];
alpha = [];
rVL = rVLEstimatorHyper(r,l1,l2,l3);        
if (isreal(rVL)&&~isnan(rVL)&&~isinf(rVL))
                        
    %% Hyper
    xtil = [x(1,:)/f;x(2,:)/f;x(3,:)];   
    chi = getXiThreePointsLinear(xtil);

    % [tanPlus tanMinus]=getTanChi(r,rVLHyper,l1,l2,l3);
    if (isreal(chi) && ~isnan(chi))
        alpha = (-xtil(3,:) + cos(chi)*sqrt(xtil(1,:).^2+xtil(2,:).^2+xtil(3,:).^2))/sin(chi);
        [~ , ~, V]=svd([xtil(1,:);-xtil(2,:);-alpha]');    
        n = V(:,end);            
        alpha = alpha*f;    
    end    
end