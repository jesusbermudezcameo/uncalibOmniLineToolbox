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

function chiEst=getXiThreePointsLinear(x_cam)

x1 = x_cam(1,1);
y1 = x_cam(2,1);
r1q = x1^2+y1^2;

x2 = x_cam(1,2);
y2 = x_cam(2,2);
r2q = x2^2+y2^2;

x3 = x_cam(1,3);
y3 = x_cam(2,3);
r3q = x3^2+y3^2;

l1 = x2*y3-x3*y2;
l2 = -x1*y3+x3*y1;
l3 = x1*y2-x2*y1;

xi = (l1+l2+l3)/(l1*sqrt(1+r1q)+l2*sqrt(1+r2q)+l3*sqrt(1+r3q));


if xi >= 1
    xi = 0.99999;
end
    
chiEst = acos(xi);



return
 
