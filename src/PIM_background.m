% Computes the PIM background model given all ChIPseq nmers
load_params;
if nearest_neighbors == false
    [~,~] = mkdir(outdirPIM);
end

%% Compute observed frequencies from all ChIPseq nmers
c1=(histc(nmersint,1:q)+pseudocount)/(size(nmersint,1)+q*pseudocount);
c2 = zeros(q,q,N,N);
for i=1:N
    for j=1:N
        c2(:,:,j,i)=reshape((histc((nmersint(:,i)-1)*q+nmersint(:,j),1:q^2)+pseudocount/q)/(size(nmersint,1)+q*pseudocount),q,q);
    end
end

%% Add correlations from independent model
J=zeros(q,q,N,N);
h=log(c1);

if nearest_neighbors == true
    % we only want to fit Ji,i+1
    tc2= false(q,q,N,N);
    for i=1:(N-1)
        tc2(:,:,i+1,i) = true;
        tc2(:,:,i,i+1) = true;
    end    
else
    % we do not want to fit Jii
    tc2= true(q,q,N,N);
    for i=1:N
        tc2(:,:,i,i) = false;
    end
end

% This is the gradient descent algorithm
EPS=1;
err1=10; % required precision in fitting h, J
while(err1>1e-5)
    
    [c1mt,c2mt,thermo]=exhaust(h,J);
    PartFunction=thermo(4);
    energytemp=energy(nmersint,h,J);
    probamodeltemp=exp(energytemp)/PartFunction;
    
    c1m=reshape(c1mt,q,N);
    c2m=reshape(c2mt,[q q N N]);
    
    % the derivatives of Dkl with respect to h and J
    % remember that d(Dkl)/d(beta) = < O >_m - < O >_r
    % where beta is the lagrangian multiplier of observale O
    tmp1=c1-c1m;
    tmp2=c2-c2m;
    
    tmp2(~tc2)=0;
    
    gg1=tmp1;
    gg2=tmp2;
    
    derivee=[gg1(:);gg2(:)];
    
    beta=0;
    dd1=0; dd2=0;
    
    dd1=gg1+beta*dd1;
    dd2=gg2+beta*dd2;
    
    % Newton's method to find minimum along the line
    g1=1;
    realEPS=EPS*.1;
    saverealEPS=realEPS;
    prevEPS=0;
    g0=sum(dd2(:).*gg2(:))+sum(dd1(:).*gg1(:));
    while(abs(g1)>1e-6)
        [c1mt,c2mt,thermo]=exhaust(h+realEPS*dd1,J+realEPS*dd2);
        PartFunction=thermo(4);
        c1m=reshape(c1mt,q,N);
        c2m=reshape(c2mt,[q q N N]);
        
        tmp1=c1-c1m;
        tmp2=c2-c2m;
        
        tmp2(~tc2)=0;
        
        gg1=tmp1;
        gg2=tmp2;
        
        g1=sum(dd2(:).*gg2(:))+sum(dd1(:).*gg1(:));
        saverealEPS=realEPS;
        realEPS=prevEPS+(realEPS-prevEPS)*g0/(g0-g1);
        if(realEPS>2000 || realEPS<0 || realEPS~=realEPS), realEPS=saverealEPS+0*g0; clear beta; break; end
        g0=g1;
        prevEPS=saverealEPS;
        fprintf('eps = %g, g0 = %g, g1 = %g\n',realEPS,g0,g1); %end
    end
    h=h+realEPS*dd1;
    J=J+realEPS*dd2;
    
    err1=sum(gg1(:).^2)+sum(gg2(:).^2);
    
    fprintf('Partfunction=%d err1=%g c1=%g c2=%g \n',PartFunction, err1,...
        sum(abs(c1(:)-c1m(:))),...
        sum(abs(c2(:)-c2m(:))))
    
end


[c1m,c2m,thermo]=exhaust(h,J);
PartFunction=thermo(4);
energytemp=energy(nmersint,h,J);
probaback=exp(energytemp)/PartFunction;

Jback = J;
hback = h;

if nearest_neighbors == false
    probabackpim = probaback;
    fprintf('Saving to %s\n',filePIMback);
    save(filePIMback,'hback','Jback','probabackpim');
else
    probabacknnm = probaback;
    fprintf('Saving to %s\n',fileNNMback);
    save(fileNNMback,'hback','Jback','probabacknnm');
end

