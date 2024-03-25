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

function dalpha_r=getDalpha_r(r,calib,modelClass)
f=calib(1);
rVL =f;

switch(modelClass)
    case 1,
        dalpha_r= 1/rVL;            
    case 2,
        chi = calib(2);
        dalpha_r = cot(chi)./sqrt(r.^2+f^2);        
    case 3,
        dalpha_r = (cot(r/f).^2-cot(r/f)*f./r+1)/f;        
    case 4,
        dalpha_r= 1/rVL;        
    case 5,
        dalpha_r = 1./sqrt(r.^2-f^2);        
    case 6,
        dalpha_r = (3*f^2-2*r.^2)./(2*(f^2-r.^2).^(3/2));        
%              dalpha_r = r.*(3*f^2-2*r.^2)./(2*(f^2-r.^2).^(3/2));           
end


return;

% Para
% ddalpha = 1/rVL;

% Hyper
% ddalpha = (f^2*cot(chi))/(f^2 + r^2)^(3/2)

% Equiangular
% ddalpha = (2*(f - r*cot(r/f))*(cot(r/f)^2 + 1))/f^2

% Orthogonal
% ddalpha = -f^2/(r^2 - f^2)^(3/2)

% Equisolid
% ddalpha = (r*(6*f^4 - 5*f^2*r^2 + 2*r^4))/(2*(f^2 - r^2)^(5/2))
