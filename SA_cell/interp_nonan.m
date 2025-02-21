function Yn=interp_nonan(Y,u,i,L,shift);
% for handling NaNs
ind=isnan(Y(:,u,i,L));
rr=1:length(Y(:,u,i,L));
noNAN=find(~ind);
Yn(:,u,i,L)=Y(:,u,i,L);
if ind(1)==1||ind(2)==1
    m=noNAN(shift);
    Yy(1:length(Y(:,u,i,L))-m)=Y(m+1:end,u,i,L);
    Yy(length(Y(:,u,i,L))-m+1:length(Y(:,u,i,L)))=Y(1:m,u,i,L);
    ind2=isnan(Yy);
    Ynew=interp1(rr(~ind2),Yy(~ind2),rr);
    Yn(m+1:end,u,i,L)=Ynew(1:length(Y(:,u,i,L))-m);
    Yn(1:m,u,i,L)=Ynew(length(Y(:,u,i,L))-m+1:length(Y(:,u,i,L)));
elseif ind(end)==1||ind(end-1)==1
    m=noNAN(shift);
    Yy(1:length(Y(:,u,i,L))-m)=Y(m+1:end,u,i,L);
    Yy(length(Y(:,u,i,L))-m+1:length(Y(:,u,i,L)))=Y(1:m,u,i,L);
    ind2=isnan(Yy);
    Ynew=interp1(rr(~ind2),Yy(~ind2),rr);
    Yn(m+1:end,u,i,L)=Ynew(1:length(Y(:,u,i,L))-m);
    Yn(1:m,u,i,L)=Ynew(length(Y(:,u,i,L))-m+1:length(Y(:,u,i,L)));
else
        Yn(:,u,i,L)=interp1(rr(~ind),Y(~ind,u,i,L),rr);
end