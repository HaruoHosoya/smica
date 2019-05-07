function [J,grad]=icaobj(params,X,nin,nout)

W0=reshape(params,nout,nin);

Wnorm=sqrt(sum(W0.^2,2));
W=bsxfun(@rdivide,W0,Wnorm);

[Y,Y1,Y2]=g(W*X);

Z=W'*Y;
U=W*Z;

J=sum(sum(Y1))+sum(sum(Z.^2))/2;

gradW=Y2*X'+Y*Z'+(Y1.*U)*X';

gradW=bsxfun(@rdivide,gradW,Wnorm)-bsxfun(@times,W,sum(W.*gradW,2)./Wnorm);

grad=gradW(:);

end

function [Y,Y1,Y2]=g(X)

c=1;
G=-sqrt(c+X.^2);
Y=X./G;
Y1=c./(G.^3);
Y2=-3*c*X./(G.^5);

end

