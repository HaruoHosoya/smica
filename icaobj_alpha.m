function [J,grad]=icaobj_alpha(params,X,nin,nout)

W0=reshape(params(1:end-nout),nout,nin);

logalpha=params(end-nout+1:end);
alpha=exp(logalpha);

Wnorm=sqrt(sum(W0.^2,2));
W=bsxfun(@rdivide,W0,Wnorm);

[Y,Y1,Y2]=g(W*X);
Ya=bsxfun(@times,Y,alpha);
Y1a=bsxfun(@times,Y1,alpha);
Y2a=bsxfun(@times,Y2,alpha);

Z=W'*Ya;
U=W*Z;

J=sum(sum(Y1a))+sum(sum(Z.^2))/2;

gradW=Y2a*X'+Ya*Z'+(Y1a.*U)*X';

gradW=bsxfun(@rdivide,gradW,Wnorm)-bsxfun(@times,W,sum(W.*gradW,2)./Wnorm);

grad_alpha=sum(Y1+Y.*U,2);
grad_logalpha=alpha.*grad_alpha;

grad=[gradW(:);grad_logalpha(:)];

end

function [Y,Y1,Y2]=g(X)

c=1;
G=-sqrt(c+X.^2);
Y=X./G;
Y1=c./(G.^3);
Y2=-3*c*X./(G.^5);

end

