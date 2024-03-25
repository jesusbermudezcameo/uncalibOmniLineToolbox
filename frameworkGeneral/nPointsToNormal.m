function n = nPointsToNormal(xHat,calib,modelClass)

r =sqrt(xHat(1,:).^2 + xHat(2,:).^2);
alpha=getAlpha(r,calib,modelClass);
sign=getSign(modelClass);

A = [xHat(1,:);sign*xHat(2,:);-alpha]';

[U S V] = svd(A);

n = V(:,end);
n = n/norm(n); 
 

return;
