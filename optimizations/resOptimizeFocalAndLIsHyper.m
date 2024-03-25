function res = resOptimizeFocalAndLIsHyper(Op,xCell,m)

f = Op(1);
chi = Op(2);
nPhi = Op(3:m+2);
nPhiVar = Op(m+3:2*m+2);

res = [];
for k=1:m
    xHat = xCell{k};
    n = [cos(nPhi(k))*cos(nPhiVar(k));cos(nPhi(k))*sin(nPhiVar(k));sin(nPhi(k))];
    
    r =sqrt(xHat(1,:).^2 + xHat(2,:).^2);
    alpha=getAlphaH(r,f,chi);    

    A = [xHat(1,:);-xHat(2,:);-alpha]';        
    res = [res;A*n];
    
end


return;
