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

function  [acceptedIndex, sizeAccepted, r_vlAcum] = uncalibratedLineExtractionDeviceDistance(xHatSel,xHatTest,gradThetaSel,nAttempts,pixelThreshold,sign,modelClass,fH)

global secondOrderFlag;

p=2; 
nPointsSel = size(xHatSel,2);
index=ceil(rand(nAttempts,p)*nPointsSel);

acceptedIndex = cell(1,size(index,1));
sizeAccepted = zeros(1,size(index,1));
r_vlAcum = zeros(1,size(index,1));

for k=1:size(index,1)
    xRand = xHatSel(:,index(k,:));     
    gradThetaRand = gradThetaSel(index(k,:));
    dxRand = [cos(gradThetaRand);sin(gradThetaRand)];
    
    if (xRand(1,1)==0&&xRand(2,1)==0)||(xRand(1,2)==0&&xRand(2,2)==0)
        disp('Error');
    end

    if secondOrderFlag == true
        switch(modelClass)
            case 1,
                r_vl = getVanishingLine2DegStereo(xRand,dxRand,sign);
            case 2,
                %% Not yet released
                
%                 if chiKnownFlag
%                     chi = chiGlobal;
%                     r_vl = getVanishingLine2DegHyper(xRand,dxRand,sign,chi);          
%                 else
%                     r_vlAux = getVanishingLine(xRand,dxRand,sign);
%                     chi = atan(r_vlAux/fH);  
%                     if chi<0.45, chi = 0.45; end;
%                     if chi>0.9, chi = 0.9; end;                
%                     r_vl = getVanishingLine2DegHyper(xRand,dxRand,sign,chi);          
%                 end
            case 3,
                r_vl = getVanishingLine2DegEq(xRand,dxRand,sign);
            case 4,
                r_vl = getVanishingLine2DegStereo(xRand,dxRand,sign);
            case 5,
                r_vl = getVanishingLine2DegOrtho(xRand,dxRand,sign);
            case 6,
                K1 = 2;
                K2 = 6;                               
                r_vl = getVanishingLine2DegModelKnown(xRand,dxRand,sign,K1,K2);
        end               
    else
        r_vl = getVanishingLine(xRand,dxRand,sign);
    end
    
    
    if secondOrderFlag == true
        r_vl = getVanishingLine2DegEq(xRand,dxRand,sign);
    else
        r_vl = getVanishingLine(xRand,dxRand,sign);
    end
    
    
    r_vlAcum(k) = r_vl;
    
    if isreal(r_vl)&&~isnan(r_vl)&&~isinf(r_vl)
                
        % 1 ParamX - 2 Hyper - 3 Equiangular - 4 StereographicX - 5 Orthogonal - 6 Equisolid   

        calib = getCalibFromRVL(r_vl,modelClass,fH);
        rTest = sqrt(xHatTest(1,:).^2+xHatTest(2,:).^2);
        iAccepted = testCandidateNoAngle(xHatSel,xHatTest,rTest,calib,modelClass,pixelThreshold);    
        sizeAccepted(k) = numel(iAccepted);
        acceptedIndex{k} = iAccepted; 
    else
        sizeAccepted(k) = 0;
        acceptedIndex{k} = [];        
    end
                                 
end




return
