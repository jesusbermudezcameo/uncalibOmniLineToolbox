function res = resOptimizeFocalAndLIsHyper(Op,xCell,m,f)

chi = Op(1);
nPhi = Op(2:m+1);
nPhiVar = Op(m+2:2*m+1);

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
