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
%     �Line extraction in uncalibrated central images with revolution symmetry�,
%     J. Bermudez-Cameo, G. Lopez-Nicolas and J. J. Guerrero,  
%     24th British Machine Vision Conference, BMVC, pp 1-11; Bristol, UK, Sept.  2013, 
% 
%     This work has been supported by the University of Zaragoza , the 
%     Spanish project VINEA DPI2012-31781, DGA-FSE(T04 and FEDER funds.
%     Jesus Bermudez-Cameo was supported by the FPU program AP2010-3849.

function r_vl = getRVLFromCalib(calib,modelClass)
% 1 ParamX - 2 Hyper - 3 Equiangular - 4 StereographicX - 5 Orthogonal - 6 Equisolid


switch(modelClass)
    case 1,               
        r_vl = calib ;    % calib = r_vl = 2*f*p
    case 2,        
        r_vl = calib(1)*tan(calib(2));
    case 3,        
        r_vl = calib*pi/2;  % calib = f = 2*r_vl/pi
    case 4,        
        r_vl = calib;   % calib = r_vl = 2*f
    case 5,        
        r_vl = calib;   %calib = f = r_vl 
    case 6,        
        r_vl = calib/sqrt(2);   %calib = f = r_vl*sqrt(2);
end