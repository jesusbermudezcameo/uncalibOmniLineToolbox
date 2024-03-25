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

function [rVlCand nArray xAccepted nNan] = avoidNanInf(r_vlArray,nArray,xAccepted)

bool1 = ~isnan(r_vlArray)&~isinf(r_vlArray);

bool2 = ~isnan(nArray(1,:))&~isinf(nArray(1,:));
bool3 = ~isnan(nArray(2,:))&~isinf(nArray(2,:));
bool4 = ~isnan(nArray(3,:))&~isinf(nArray(3,:));
iNC = find(bool1 & bool2 & bool3 & bool4);

rVlCand = r_vlArray(iNC);
nArray = nArray(:,iNC);
xAccepted = xAccepted(iNC);

nNan = numel(r_vlArray)-numel(iNC);

