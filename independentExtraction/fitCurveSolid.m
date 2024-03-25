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

function [n,rVL,alpha]=fitCurveSolid(x,r,l1,l2,l3)

rVL = getrVLSolidPoly(l1,l2,l3,r(1),r(2),r(3));
f = rVL*sqrt(2);

if numel(rVL) ~= 1
    errorFlag = 1;       
    rVLEst = rVLEstimatorSolid(r,l1,l2,l3);   
    if numel(rVL)>1
        absErr = abs(rVL-rVLEst);
        rVL = rVL(absErr==min(absErr));
        f = rVL*sqrt(2);
    else
%         disp('Error in polynomial estimation - fitCurveSolid.m/43');  
        rVL = rVLEst; f = rVL*sqrt(2);
    end            
end


if ((numel(rVL)==1)&&isreal(rVL)&&~isnan(rVL)&&~isinf(rVL))   
    alpha = (2*r.^2-f^2)./sqrt(f^2-r.^2)/2;            
    [~,~,V]=svd([x(1,:);x(2,:);-alpha]');    
    n = V(:,end);         
else    
    n = [];
    alpha = [];
    if (numel(rVL)~=1), rVL = 1i;end; 
end



return;
% n = [];
% alpha = [];
% 
% rVL = rVLEstimatorSolid(r,l1,l2,l3);
% if (isreal(rVL)&&~isnan(rVL)&&~isinf(rVL))
%     Op = rVL*sqrt(2);
%     options = optimset(optimset('lsqnonlin'),'MaxIter',5000,'MaxFunEvals',20000,'TolX',0.00000001,'Display','off');
%     f=lsqnonlin('resSolidRVL',Op,[],[],options,l1,l2,l3,r);
%     rVL = f/sqrt(2);  
%     
%     alpha = (2*r.^2-f^2)./sqrt(f^2-r.^2)/2;            
%     [U S V]=svd([x(1,:);x(2,:);-alpha]');    
%     n = V(:,end);                 
% end         
% 
% return;
