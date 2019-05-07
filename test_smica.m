function [A,W,alpha]=test_smica(data,nout)

% dataset32=h5read('~/datasets/imagenet/imagenet10k/processed2/patch32_100k/patches000001.h5','/data');
% dataset16=dataset32(9:24,9:24,:);
% dataset16=reshape(dataset16,16*16,size(dataset16,3));
% data=bsxfun(@minus,dataset16,mean(dataset16,1));

% nin=size(data,1)-5;
nin=size(data,1);
% nout=40;

[E,~,D]=pca(data','NumComponents',nin);
fprintf('Eigen value max=%f min=%f\n',D(1),D(nin));

V=E*diag(D(1:nin).^(-1/2));
Vi=E*diag(D(1:nin).^(1/2));
data=V'*data;

U0=randn(nout,nin);
beta0=randn(nout,1)*0.01;
% params0=[U0(:); beta0(:)];
params0=[U0(:)];


objopts.fix_alpha=true;

obj=@(params)icaobj(params,data,nin,nout);

%% Using fminunc

opts=optimoptions('fminunc','Algorithm','trust-region','GradObj','on',...
    'DerivativeCheck','on','Display','iter',...
    'TolFun',1e-3,'TolX',1e-3,'MaxFunEvals',1000);

tic
[params_opt,Jopt,eval]=fminunc(obj,params0(:),opts);
toc

%% Using minFunc

% opts.Method='lbfgs';
% opts.maxIter=1000;	    % Maximum number of iterations of L-BFGS to run 
% opts.display='iter';
% 
% tic
% [params_opt,Jopt] = minFunc(obj,params0(:),opts);
% toc

%%

U=reshape(params_opt(1:nout*nin),nout,nin);
alpha=exp(params_opt(nout*nin+1:end));
U=bsxfun(@rdivide,U,sqrt(sum(U.^2,2)));

W=U*V';
A=Vi*U';


end
