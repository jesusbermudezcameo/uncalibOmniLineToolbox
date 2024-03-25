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


function [r_vl, n, v] = getVanishingLine2DegModelKnown(xHat,dx,s,K1,K2)

r = sqrt(xHat(1,:).^2+xHat(2,:).^2);
tau = dx(2,:).*xHat(1,:)-dx(1,:).*xHat(2,:);

% N=[[xHat(1,:); s*xHat(2,:);-r*K1; ones(1,numel(r))*K1] [-dx(2,:); s*dx(1,:); K1*tau./r; zeros(1,numel(r))]]';
N = [[xHat(1,:); s*xHat(2,:);-r*(K1-K2); ones(1,numel(r))*(K1-K2/2);-K2/2*r.^2]    ,[-dx(2,:); s*dx(1,:); (K1-K2)*tau./r; zeros(1,numel(r));K2*tau]]';
[~,~,V] = svd(N);
v = V(:,end);

r_vl = v(3)/v(5);
n = v(1:3)/norm(v(1:3));

return