% Computes the PIM (Pairwise Interaction Model) given an initial PWM and a
% list of ChIPseq nmers. Change the 'nearest_neighbors' variable to true to
% compute NNM.
load_params;
if nearest_neighbors == false
    [~,~] = mkdir(outdirPIM);
end

%% Computing background model for scoring
if nearest_neighbors == false
    if (~exist(filePIMback,'file'))
        fprintf('Computing PIM background model...\n')
        PIM_background;
    else
        fprintf('Loading PIM background model...\n')
        load(filePIMback, 'probabackpim')
    end
    probaback = probabackpim;
    fileout = filePIM;
else
    fprintf('Restricting PIM to nearest neighbors...\n')
    if (~exist(fileNNMback,'file'))
        fprintf('Computing NNM background model...\n')
        NNM_background;
    else
        fprintf('Loading NNM background model...\n')
        load(fileNNMback, 'probabacknnm')
    end
    probaback = probabacknnm;
    fileout = fileNNM;
end



%% Computing PWM (if not already done)
% used ot get an initial set of sequences and initial h, J
if (~exist(filePWM,'file'))
    fprintf('Computing PWM model...\n')
    PWM_model
else
    fprintf('Loading PWM model...\n')
    load(filePWM)
end

%% Adding correlations to the PWM
fprintf('Introducing correlations...\n');

% we keep track of all quantities
SaveJ=cell(1);
Saveh=cell(1);
SavePartFunction=cell(1);
Savealloc1=cell(1);
Savealloc2=cell(1);
SaveDKL=cell(1);
Savescore=cell(1);
SaveBIC=cell(1);
Savediffsites=zeros(1);
Count=1;


% Now we begin the main loop: we constrain the model from a current set of
% sites and predict a new set of sites with it. The loop stops when the
% model predicts the same sites than the previous model or when it has a
% worse BIC
prevBIC=Inf;
while(1)
    fprintf('Iteration # %d\n',Count);
    
    prevseqs=seqs;
    seqsint=double(nt2int(seqs));
    M=size(seqsint,1);
    q=4;
    
    % Building the PWM model (c1)
    c1=(histc(seqsint,1:q)+pseudocount)/(M+q*pseudocount);
    for i=1:N
        for j=1:N
            c2(:,:,j,i)=reshape((histc((seqsint(:,i)-1)*q+seqsint(:,j),1:q^2)+pseudocount/q)/(M+q*pseudocount),q,q);
            if(i==j)
                c2i(:,:,j,i)=diag(c1(:,i));
                c2(:,:,j,i)=diag(c1(:,i));
            else
                c2i(:,:,j,i)=c1(:,j)*c1(:,i)';
            end
        end
    end
    
    % No interactions (PWM)
    hinde=log(c1);
    Jinde=zeros(q,q,N,N);
    
    
    [useqs,~,c]=unique(num2cell(seqs,2));
    useqs=cell2mat(useqs);
    useqsint=double(nt2int(useqs));
    useqslabel=sum((useqsint-1).*repmat(q.^(0:N-1),size(useqsint,1),1),2)+1;
    
    energyinde=energy(useqsint,hinde,Jinde);
    rawproba=histc(c,1:max(c))/numel(c);
    
    tochange=false(q,N);
    for i=1:N
        tmp=find(c1(:,i)>0);
        tochange(tmp,i)=true;
    end
    
    tochange2=false(q,q,N,N);
    for i=1:N
        for j=1:N
            tochange2(:,:,j,i)=logical(double(tochange(:,j))*double(tochange(:,i)'));
        end
    end
    
    PartFunction=1;
    
    c1m=c1;
    c2m=c2i;
    
    alloc1= false(N,1);
    alloc2= false(N,N);
    
    tc1= false(q,N);
    tc2= false(q,q,N,N);
    
    h=hinde;
    J=Jinde;
    
    EPS=1;
    
    saveJ=cell(1);
    saveh=cell(1);
    savePartFunction=cell(1);
    savealloc1=cell(1);
    savealloc2=cell(1);
    saveDKL=zeros(1);
    savescore=zeros(1);
    saveBIC=zeros(1);
    count=1;
    
    % Fitting the parameters by minimizing dkl
    oldloc=0;
    while(1)
        energytemp=energy(useqsint,h,J);
        probamodeltemp=exp(energytemp)/PartFunction;
        saveJ{count}=J;
        saveh{count}=h;
        savePartFunction{count}=PartFunction;
        savealloc2{count}=alloc2;
        savealloc1{count}=alloc1;
        savetc2{count}=tc2;
        savetc1{count}=tc1;
        saveDKL(count)=sum(rawproba.*log2(rawproba./probamodeltemp));
        savescore(count)=saveDKL(count)+(3*N+sum(tc2(:))/2)/2*log2(M)/M;
        saveBIC(count) = 2 * M * saveDKL(count) + (3 * N + sum(tc2(:)) / 2) * log2(M);
        fprintf('%d interactions added, DKL %g BIC %g\n',count,saveDKL(count),saveBIC(count));
        
        % "pvalue" on two point correlations to select the one to be fitted
        % c2's are the observed two-points counts
        % c2m's are the modeled two-points probabilities
        ppdf=binopdf(floor(c2*M),M+0*c2,c2m)./binopdf(floor(c2*M),M+0*c2,c2);
        for i=1:N
            ppdf(:,:,i,i)=1;
        end
        ppdf(abs(c2-c2m)<=2/M & abs(c2)<10/M)=1;
        for i=1:N
            for j=1:N
                if(sum(sum(tc2(:,:,i,j)))>=9)
                    ppdf(:,:,i,j)=1;
                end
                for a=1:q
                    for b=1:q
                        if tc2(a,b,i,j) == 1
                            ppdf(a,b,i,j) = 1; % do not learn old interactions again
                        end
                    end
                end
            end
        end
        
        % this piece is to keep only nearest neighbors
        if nearest_neighbors == true
            for i=1:N
                for j=1:N
                    if abs(j-i) > 1
                        ppdf(:,:,i,j)=1;
                    end
                end
            end
        end
        
        [m,loc]=min(ppdf(:));
        % criterium to stop the addition of parameters: the BIC increases
        % on average during the last averageextent iterations. This is only
        % for illustration purposes, to go faster one should use
        % averageextent = 1;
        averageextent=20;
        meanderiv=0;
        if (count>averageextent)
            meanderiv=mean(diff(savescore(count-averageextent:count)));
            if (meanderiv>0)
                [a,b,i,j]=ind2sub(size(ppdf),loc);
                %                 fprintf('i=%d j=%d a=%d b=%d log(m)=%g meanderiv=%g\n\n',i,j,a,b,log(m),meanderiv);
                break;
            end
        end
        
        [a,b,i,j]=ind2sub(size(ppdf),loc);
        alloc2(i,j)=true; alloc2(j,i)=true;
        alloc1(i)=true; alloc1(j)=true;
        tc2(a,b,i,j)=true; tc2(b,a,j,i)=true;
        tc1(a,i)=true; tc1(b,j)=true;
        
        if(loc~=oldloc)
            clear beta;
            count=count+1;
        end
        %         fprintf('i=%d j=%d a=%d b=%d log(m)=%g meanderiv=%g\n\n',i,j,a,b,log(m),meanderiv);
        
        oldloc=loc;
        
        tmpDKL=10;
        % This is the gradient descent algorithm
        err1=10; % precision in fitting h, J
        while(err1>1e-5)
            [c1mt,c2mt,thermo]=exhaust(h(:,alloc1),J(:,:,alloc1,alloc1));
            PartFunction=thermo(4);
            energytemp=energy(useqsint,h,J);
            probamodeltemp=exp(energytemp)/PartFunction;
            
            oldDKL=tmpDKL;
            
            tmpDKL=sum(rawproba.*log2(rawproba./probamodeltemp));
            c1m(:,alloc1)=reshape(c1mt,q,length(find(alloc1)));
            c2m(:,:,alloc1,alloc1)=reshape(c2mt,[q q length(find(alloc1)) length(find(alloc1))]);
            
            tmp1=c1-c1m;
            tmp2=c2-c2m;
            
            tmp1(~(tc1 & tochange))=0;
            tmp2(~(tc2 & tochange2))=0;
            
            gg1=tmp1(:,alloc1);
            gg2=tmp2(:,:,alloc1,alloc1);
            
            derivee=[gg1(:);gg2(:)];
            
            erreur=sum(derivee.^2);
            if(exist('beta','var') && length(find(alloc1))>20) % we avoid convergence problems
                newscalar=sum(derivee.^2);
                beta=max(0,(newscalar-sum(oldgg2(:).*gg2(:))-sum(oldgg1(:).*gg1(:)))/oldscalar);
                %beta=newscalar/oldscalar;
                oldscalar=newscalar;
            else
                beta=0;
                dd1=0; dd2=0;
                oldscalar=sum(derivee.^2);
            end
            oldgg2=gg2;oldgg1=gg1;
            
            dd1=gg1+beta*dd1;
            dd2=gg2+beta*dd2;
            
            g1=1;
            realEPS=EPS;
            %while(abs(log(prevEPS)-log(realEPS))>1e-10)
            if(length(find(alloc1))>8)
                realEPS=EPS*.1;
                saverealEPS=EPS*.1;
                prevEPS=0;
                g0=sum(dd2(:).*gg2(:))+sum(dd1(:).*gg1(:));
                % Newton's method to find minimum along the line
                while(abs(g1)>1e-6)
                    [c1mt,c2mt,thermo]=exhaust(h(:,alloc1)+realEPS*dd1,J(:,:,alloc1,alloc1)+realEPS*dd2);
                    PartFunction=thermo(4);
                    c1m(:,alloc1)=reshape(c1mt,q,length(find(alloc1)));
                    c2m(:,:,alloc1,alloc1)=reshape(c2mt,[q q length(find(alloc1)) length(find(alloc1))]);
                    
                    tmp1=c1-c1m;
                    tmp2=c2-c2m;
                    
                    tmp1(~(tc1 & tochange))=0;
                    tmp2(~(tc2 & tochange2))=0;
                    
                    gg1=tmp1(:,alloc1);
                    gg2=tmp2(:,:,alloc1,alloc1);
                    
                    g1=sum(dd2(:).*gg2(:))+sum(dd1(:).*gg1(:));
                    saverealEPS=realEPS;
                    realEPS=prevEPS+(realEPS-prevEPS)*g0/(g0-g1);
                    if(realEPS>2000 || realEPS<0 || realEPS~=realEPS), realEPS=saverealEPS+0*g0; clear beta; break; end
                    %if(length(find(alloc1))>8),
                    %fprintf('eps = %g, g0 = %g, g1 = %g\n',realEPS,g0,g1); %end
                    g0=g1;
                    prevEPS=saverealEPS;
                end
            end
            h(:,alloc1)=h(:,alloc1)+realEPS*dd1;
            J(:,:,alloc1,alloc1)=J(:,:,alloc1,alloc1)+realEPS*dd2;
            
            err1=sum(gg1(:).^2)+sum(gg2(:).^2);
            
            %err1=oldDKL-tmpDKL;
            
            %             fprintf('Partfunction=%d err1=%g erreur=%g c1=%g c2=%g alloc1=%d alloc2=%d tc2=%d tmpDKL=%g score=%g\n',PartFunction, err1,erreur,...
            %                sum(abs(c1(tochange)-c1m(tochange))),...
            %                sum(abs(c2(tochange2)-c2m(tochange2))),...
            %                sum(alloc1),sum(sum(tril(alloc2,-1))),...
            %                sum(tc2(:))/2,tmpDKL,tmpDKL+(3*N+sum(tc2(:))/2)/2*log2(M)/M);
        end
    end
    
    if min(saveBIC) > prevBIC
        fprintf('New model is worst than previous iteration, stopping.\n');
        break
    end
    prevBIC = min(saveBIC);
    
    [~,ll]=min(savescore);
    J=saveJ{ll};
    h=saveh{ll};
    PartFunction=savePartFunction(ll);
    
    SaveJ{Count}=saveJ;
    Saveh{Count}=saveh;
    SavePartFunction{Count}=savePartFunction;
    Savealloc2{Count}=savealloc2;
    Savealloc1{Count}=savealloc1;
    Savetc2{Count}=savetc2;
    Savetc1{Count}=savetc1;
    SaveDKL{Count}=saveDKL;
    Savescore{Count}=savescore;
    SaveBIC{Count}=saveBIC;
    
    
    % Now we find nmers probas to find the sequences predicted by the refined model
    PartFunction=thermo(4);
    energytemp = energy(nmersint,h,J);
    probamodeltemp = exp(energytemp) / PartFunction;
    
    scoresites = log(probamodeltemp ./ probaback);
    
    [sproba,i]=sort(scoresites,'descend');
    for nmax=0:length(i)
        if length(unique(enhs(i(1:nmax)))) > Nhalf
            Nmax=nmax;
            break;
        end
    end
    
    seqs=nmers(i(1:Nmax),:);
    
    diffsites = size(setdiff(seqs, prevseqs, 'rows'),1);
    fprintf('%d site(s) differ with previous set of sequences\n', diffsites);
    Savediffsites(Count)=diffsites;
    
    Count=Count+1;
    
    fprintf('Saving to %s\n',fileout);
    save(fileout, 'h', 'J', 'PartFunction', 'seqs', 'Save*');
    
    
    if diffsites==0 && length(seqs) == length(prevseqs)
        fprintf('Converged!\n');
        break
    end
end


fprintf('Exit normally!\n');
