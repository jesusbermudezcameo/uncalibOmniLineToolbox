function  Op_optim = optimFocalAndLIsHyper(xAcceptedAcum,nAcum,f,chi)

binSel = sum(imag(nAcum)==0)==3;
nAcum = nAcum(:,binSel);
xAcceptedAcum = xAcceptedAcum(binSel);

if numel(xAcceptedAcum)==0
    Op_optim = chi;
    return;
end


nPhi = asin(nAcum(3,:));
nPhiVar = atan2(nAcum(2,:),nAcum(1,:));

Op = [chi nPhi nPhiVar];
limInf = [0 -pi/2*ones(1,numel(xAcceptedAcum)) -pi*ones(1,numel(xAcceptedAcum))];
limSup = [pi/2 pi/2*ones(1,numel(xAcceptedAcum)) pi*ones(1,numel(xAcceptedAcum))];
options = optimset(optimset('lsqnonlin'),'MaxIter',5000,'MaxFunEvals',20000,'TolX',0.00000001,'Display','off');      
Op_optim = lsqnonlin('resOptimizeFocalAndLIsHyper',Op,limInf,limSup,options,xAcceptedAcum,numel(xAcceptedAcum),f);    


return;