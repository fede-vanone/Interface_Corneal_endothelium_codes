function [w,x]=clenshaw_curtis_p2(N,p)
% like clenshaw_curtis_p2.m but doesn't need the PDE toolbox.


jW=0:1:N-1; jW=jW(:);

%integration weights p=1

IWp1=diag(INTweights(N,2));



% I type moments

TM=cos(repmat(jW,[1,N]).*repmat(0:N-1,[N,1])*(pi/(N-1)));

pM=repmat(p,[1,N]);

IWp1M=repmat(IWp1,[1,N]);

momsI=sum(TM.*pM.*IWp1M).';

% II type moments 

momsII=0.5*(momsI(1:end-2)-momsI(3:end));



% integration weights p=p

theta=(1:N-2)'*pi/(N-1); xx=cos(theta);

w=((2*sin(theta)/(N-1)).*dst2(momsII))./(1-xx.^2);

w1=(2*sum(momsI)-momsI(1)-momsI(end))/(2*(N-1));

wn=momsI(1)-w1-sum(w);

x=[1;xx;-1]; w=[w1;w;wn];


end


function IW=INTweights(N,boxsize)



nW=0:1:N-1;

jW=0:1:N-1;



bW=ones(1,N);

bW(1)=0.5;

bW(N)=0.5;

cW=2*bW;

bW=bW/(N-1);

S=cos(nW(3:N)'*jW*(pi/(N-1)));

IW=boxsize/2*diag(bW.*[(2+(cW(3:N).*((1+(-1).^nW(3:N))./(1-nW(3:N).^2)))*S)]);

end



function b=dst2(a,n)



narginchk(1,2);



if min(size(a))==1

    if size(a,2)>1

        do_trans = 1;

    else

        do_trans = 0;

    end

    a = a(:);

else

    do_trans = 0;

end

if nargin==1

  n = size(a,1);

end

m = size(a,2);



% Pad or truncate a if necessary

if size(a,1)<n

  aa = zeros(n,m);

  aa(1:size(a,1),:) = a;

else

  aa = a(1:n,:);

end



y=zeros(2*(n+1),m);

y(2:n+1,:)=aa;

y(n+3:2*(n+1),:)=-flipud(aa);

yy=fft(y);

b=yy(2:n+1,:)/(-2*sqrt(-1));



if isreal(a), b = real(b); end

if do_trans, b = b.'; end



% LocalWords:  Nordmark

end

