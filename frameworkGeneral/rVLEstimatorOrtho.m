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

function rVLOrtho = rVLEstimatorOrtho(r,l1,l2,l3)

s1 = l1+l2+l3;
s5 = l1*r(1)^2+l2*r(2)^2+l3*r(3)^2;


%% Orthogonal
rVLOrtho =  sqrt(s5/s1/2);

a = s1;
b = -0.5*s5;
c = -1/8*(l1*r(1)^4+l2*r(2)^4+l3*r(3)^4);

rVLOrtho2 = sqrt((-b+sqrt(b^2-4*a*c))/2/a);
rVLOrtho3 = sqrt((-b-sqrt(b^2-4*a*c))/2/a);

if isreal(rVLOrtho2) && ~isreal(rVLOrtho3) 
    rVLOrtho = rVLOrtho2;    
else
    if isreal(rVLOrtho3) && ~isreal(rVLOrtho2)
        rVLOrtho = rVLOrtho3;                       
    end    
end
return;