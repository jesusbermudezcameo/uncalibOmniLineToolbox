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

function [n, rVL, alpha, l1, l2, l3, r] = fitCurve(x,modelClass,fH)


[l1, l2, l3, r] = getrl(x);
% 1 ParamX - 2 Hyper - 3 Equiangular - 4 StereographicX - 5 Orthogonal - 6 Equisolid
switch(modelClass)
    case 1,
        [n, rVL, alpha]=fitCurvePara(x,r,l1,l2,l3);       
    case 2,
        [n, rVL, alpha]=fitCurveHyper(x,r,l1,l2,l3,fH);     
    case 3,
        [n, rVL, alpha]=fitCurveEqui(x,r,l1,l2,l3);    
    case 4,
        [n, rVL, alpha]=fitCurveStereo(x,r,l1,l2,l3);                                
    case 5,
        [n, rVL, alpha]=fitCurveOrtho(x,r,l1,l2,l3);  
    case 6,
        [n, rVL, alpha]=fitCurveSolid(x,r,l1,l2,l3);
end



return;