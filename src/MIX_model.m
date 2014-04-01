% Computes the PWM mixture model given an initial PWM and a list of ChIPseq
% nmers
load_params;
[~,~] = mkdir(outdirMIX);

%% Computing PWM (if not already done)
% used ot get an initial set of sequences and initial h, J
if (~exist(filePWM,'file'))
    fprintf('Computing PWM model...\n')
    PWM_model
else
    fprintf('Loading PWM background model...\n')
    load(filePWMback)
    fprintf('Loading PWM model...\n')
    load(filePWM)
end

%% Use K-MEANS to define sets of related sites
fprintf('K-MEANS clustering...\n');
Count=1;
Savediffsites=zeros(1);
Saveh=cell(1);
Saveseqs=cell(1);
Saveassignments=cell(1);
Savecenters=cell(1);
Saveclusprob=cell(1);
SaveDKL=cell(1);
Savescore=cell(1);
SaveBIC=cell(1);
Saveprobamodel=cell(1);

prevBIC=Inf;

% Now we begin the main loop: we constrain the model from a current set of
% sites and predict a new set of sites with it. The loop stops when the
% model predicts the same sites than the previous model or when it has a
% worse BIC
while(1)

%    fprintf('Iteration # %d\n',Count);
   
   prevseqs = seqs;
   seqsint = double(nt2int(seqs));
   seqsbin = nt2bin(seqs); % transforms to 4*L binary sequence
   N=size(seqsint,2);
   M=size(seqsint,1);
   q=4;
    
   % Building the PWM model (c1)
   K=1;
   c1=(histc(seqsint,1:q)+pseudocount)/(M+q*pseudocount); 
   hinde=log(c1);

   [useqs,~,c]=unique(num2cell(seqs,2));
   useqs=cell2mat(useqs);
   useqsint=double(nt2int(useqs));
    
   energyinde=energy(useqsint,hinde,Jinde);
   probamodeltemp=exp(energyinde);
   rawproba=histc(c,1:max(c))/numel(c);
    
   h=cell(1);
   h{1}=hinde;
   J=Jinde;

   saveh=cell(1);
   saveDKL=zeros(1);
   savescore=zeros(1);
   saveBIC=zeros(1);
   saveIDX=cell(1);
   saveC=cell(1);
   saveclusprob=cell(1);
   saveprobamodel=cell(1);
   count=1;

   saveprobamodel{count}=probamodeltemp;
   saveclusprob{count}=1;
   saveIDX{count}=ones(M,1);
   saveC{count}=mean(seqsint,1);
   saveh{count}=h;
   saveDKL(count)=sum(rawproba.*log2(rawproba./probamodeltemp));
   savescore(count)=saveDKL(count)+3*N/2*log2(M)/M;
   saveBIC(count) = 2 * M * saveDKL(count) + 3 * N * log2(M);
   fprintf('K= %d DKL %g BIC %g\n',K,saveDKL(count),saveBIC(count));

   count=count+1;

   % Building mixture model
   Nrep = 20; % number of random initial cluster centers. Increase for precision.
   for (K=2:20)

      warning off;
      [IDX,C]=kmeans(seqsbin,K,'Replicates',Nrep,'Emptyaction','singleton','Distance','Hamming');
      warning on;
      saveIDX{count}=IDX;
      saveC{count}=C;

      probanmers=zeros(1);
      probamodel=zeros(1);
      c1=zeros(4,N);
      h=cell(1);
      for (k=1:K)
            clusseqsint=seqsint(IDX==k,:);
            clusM=size(clusseqsint,1);
            clusprob=clusM/M;
            saveclusprob{count}(k)=clusprob;
            % PWM
            c1=(histc(clusseqsint,1:q)+pseudocount)/(clusM+q*pseudocount);
            h{k}=log(c1);
            energymodel=energy(useqsint,h{k},Jinde);
            probamodel=probamodel+clusprob*exp(energymodel);
      end

      saveprobamodel{count}=probamodel;
      
      saveh{count}=h;
      saveDKL(count)=sum(rawproba.*log2(rawproba./probamodel));
      savescore(count)=saveDKL(count)+K*N*3/2*log2(M)/M;
      saveBIC(count) = 2 * M * saveDKL(count) + 3 * N * K * log2(M);
      fprintf('K= %d DKL %g BIC %g\n',K,saveDKL(count),saveBIC(count));

      count=count+1;

   end

   if min(saveBIC) > prevBIC
      fprintf('New model is worst than previous iteration, stopping.\n');
      break
  end
  prevBIC = min(saveBIC);

  [~,ll]=min(savescore);
  probamodel=saveprobamodel{ll};
  Saveprobamodel{Count}=saveprobamodel{ll};
  Saveassignments{Count}=saveIDX{ll};
  Savecenters{Count}=saveC{ll};
  Saveclusprob{Count}=saveclusprob{ll};
  
  % the model reads sum_k clusprob_k * exp(h_k)
  h=saveh{ll};
  clusprob = saveclusprob{ll}; % clusters probabilities

  Saveh{Count}=saveh;
  SaveDKL{Count}=saveDKL;
  Savescore{Count}=savescore;
  SaveBIC{Count}=saveBIC;


  probanmers=zeros(1);
  for (k=1:ll)
         energynmers=energy(nmersint,h{k},Jinde);
         probanmers=probanmers+saveclusprob{ll}(k)*exp(energynmers);
   end
   scoresites = log(probanmers) - energybackinde;

   [sproba,i]=sort(scoresites,'descend');
   for nmax=1:length(i)
      if length(unique(enhs(i(1:nmax)))) > Nhalf
          Nmax=nmax;
          break;
      end
   end

   seqs=nmers(i(1:Nmax),:);

   diffsites = size(setdiff(seqs, prevseqs, 'rows'),1);
   fprintf('Iteration %d: K=%d, %d site(s) differ\n', Count, ll, diffsites);

   Saveseqs{Count}=seqs;
   Savediffsites(Count)=diffsites;

   Count=Count+1;

   save(fileMIX,'Save*','h','clusprob','seqs');


   if diffsites==0 && length(seqs) == length(prevseqs)
      fprintf('Converged!\n');
      break
  end

end

fprintf('Exit normally!\n');
