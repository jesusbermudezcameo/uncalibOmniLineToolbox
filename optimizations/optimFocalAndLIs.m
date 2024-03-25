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

function  Op_optim = optimFocalAndLIs(xAcceptedAcum,nAcum,f,modelClass)

sign = getSign(modelClass);

binSel = sum(imag(nAcum)==0)==3;
nAcum = nAcum(:,binSel);
xAcceptedAcum = xAcceptedAcum(binSel);

if numel(xAcceptedAcum) == 0
    Op_optim = f;
    return;
end

nPhi = asin(nAcum(3,:));
nPhiVar = atan2(nAcum(2,:),nAcum(1,:));

Op = [f nPhi nPhiVar];
options = optimset(optimset('lsqnonlin'),'MaxIter',5000,'MaxFunEvals',20000,'TolX',0.00000001,'Display','off');     
leftLimit = [0 -pi/2*ones(1,numel(nPhi)) -pi*ones(1,numel(nPhi))];
rightLimit = [10000 pi/2*ones(1,numel(nPhi)) pi*ones(1,numel(nPhi))];
Op_optim = lsqnonlin('resOptimizeFocalAndLIs',Op,[],[],options,xAcceptedAcum,numel(xAcceptedAcum),modelClass,sign);    

if isreal(Op) && ~isreal(Op_optim)
    Op_optim = Op;
end
    

return;