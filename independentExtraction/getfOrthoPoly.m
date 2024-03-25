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

function fOrtho =  getfOrthoPoly(l1,l2,l3,r1,r2,r3)

% terms_f = [ f^4, f^2, 1];

c = [
                                                                                            - l1^4 + 2*l1^2*l2^2 + 2*l1^2*l3^2 - l2^4 + 2*l2^2*l3^2 - l3^4;
 2*l1^4*r1^2 - 2*l1^2*l2^2*r1^2 - 2*l1^2*l2^2*r2^2 - 2*l1^2*l3^2*r1^2 - 2*l1^2*l3^2*r3^2 + 2*l2^4*r2^2 - 2*l2^2*l3^2*r2^2 - 2*l2^2*l3^2*r3^2 + 2*l3^4*r3^2;
                                               - l1^4*r1^4 + 2*l1^2*l2^2*r1^2*r2^2 + 2*l1^2*l3^2*r1^2*r3^2 - l2^4*r2^4 + 2*l2^2*l3^2*r2^2*r3^2 - l3^4*r3^4];




f2A = (-c(2)+sqrt(c(2)^2-4*c(1)*c(3)))/2/c(1);
f2B = (-c(2)-sqrt(c(2)^2-4*c(1)*c(3)))/2/c(1);

f1 = sqrt(f2A);
f2 = -sqrt(f2A);
f3 = sqrt(f2B);
f4 = -sqrt(f2B);

fOrtho = [f1 f2 f3 f4];
r = [r1;r2;r3];
res(1) = resOrthoRVL(f1,l1,l2,l3,r); solReal(1) = isreal(res(1));
res(2) = resOrthoRVL(f2,l1,l2,l3,r); solReal(2) = isreal(res(2));
res(3) = resOrthoRVL(f3,l1,l2,l3,r); solReal(3) = isreal(res(3));
res(4) = resOrthoRVL(f4,l1,l2,l3,r); solReal(4) = isreal(res(4));

binZeros = abs(res)<0.01;

solPositive = fOrtho>0;

fOrtho = fOrtho(solPositive & solReal & binZeros);
return;



% A = (l1 + l2 - l3)*(l1 - l2 + l3)*(l2 - l1 + l3)*(l1 + l2 + l3)
% B/2 = (l1^4 -l1^2*(l2^2 + l3^2))*r1^2 + (l2^4 -l2^2*(l1^2 + l3^2))*r2^2 + (l3^4-l3^2*(l1^2 + l2^2) )*r3^2
% C = (l1*r1 + l2*r2 + l3*r3)*(l1*r1 + l2*r2 - l3*r3)*(l1*r1 - l2*r2 + l3*r3)*(l2*r2 - l1*r1 + l3*r3)



