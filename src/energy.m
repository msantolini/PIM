function [res,hcontr,Jcontr,hfelt]=energy(seqs,h,J)
res=0;
N=size(h,2);
q=size(h,1);
if(nargout>1 && size(seqs,1)>1)
    error('only one sequence if one wants energy contribs and felt field');
end
hcontr=zeros(1,N);
for i=1:N
    res=res+h(seqs(:,i),i);
    if(nargout>1)
        hcontr(i)=h(seqs(:,i),i);
    end
end
J=reshape(J,[q^2,N,N]);

if(nargout>1)
    hfelt=h;
end


for i=2:N
    for j=1:i-1
        res=res+J((seqs(:,i)-1)*q+seqs(:,j),j,i);
        if(nargout>1)
            Jcontr(i,j)=J((seqs(:,i)-1)*q+seqs(:,j),j,i);
            Jcontr(j,i)=Jcontr(i,j);
            hfelt(:,i)=hfelt(:,i)+J((0:3)'*q+seqs(:,j),j,i);
            hfelt(:,j)=hfelt(:,j)+J((seqs(:,i)-1)*q+(1:4)',j,i);
        end
    end
end

if(nargout>1)
    for i=1:N
        hfelt(:,i)=hfelt(:,i)-hfelt(seqs(:,i),i);
    end
end
end