% Computes the PWM model given an initial PWM and a list of ChIPseq nmers
load_params;
[~,~] = mkdir(outdirPWM);


%% Computing background PWM
if (~exist(filePWMback,'file'))
    fprintf('Computing PWM background model...\n')
    PWM_background;
else
    fprintf('Loading PWM background model...\n')
    load(filePWMback)
end

%% Changing initial PWM size to N=12
c1 = c1init; 
if (lenc1 ~= N)
   fprintf('Changing initial PWM to size %d\n',N);
   infopos = sum(c1 .* log(c1 ./ repmat(meanc1backinde, 1, lenc1)), 1);
   meanpos = 0;
   for i=1:lenc1
      meanpos = meanpos + i * infopos(i) / sum(infopos);
   end
   meanpos = floor(meanpos);
   initpos = meanpos - floor((N - 1) / 2);
   finpos = initpos + N - 1;

   finc1 = [];
   for i = initpos:finpos
      if (i < 1 || i > lenc1)
         finc1 = [finc1 meanc1backinde]; 
      else
         finc1 = [finc1 c1(:, i)]; 
      end
   end
   c1 = finc1;
end

%% Now we find sites on ChIP peaks for this initial PWM
Jinde=zeros(q,q,N,N);
hindeinit=log(c1);
energymodel = energy(nmersint,hindeinit,Jinde);
scoresites = energymodel - energybackinde;

[sproba,i]=sort(scoresites,'descend');
for nmax=1:length(i)
      if length(unique(enhs(i(1:nmax)))) > Nhalf
         Nmax=nmax;
         break;
  end
end
seqs=nmers(i(1:Nmax),:);
seqsint=double(nt2int(seqs));
M = size(seqsint,1);

initseqs = seqs;

%% Iterative procedure to refine the initial PWM given these sites
fprintf('PWM centering and refinement...\n');
Count=1;
Savediffsitespwm=zeros(1);
Savehindepwm=cell(1);
Saveseqsinde=cell(1);
prevmeanpos = floor(N / 2);
while(1)
   prevseqs=seqs;
   seqsint=double(nt2int(seqs));
   M = size(seqsint,1);

   % Building the PWM model (c1)
   c1=(histc(seqsint,1:q)+pseudocount)/(M+q*pseudocount); 

   % Centering PWM
   infopos = sum(c1 .* log(c1 ./ c1back), 1);
   meanpos = 0;
   for i=1:N 
      meanpos = meanpos + i * infopos(i) / sum(infopos);
   end
   meanpos = floor(meanpos);

   if (meanpos ~= prevmeanpos)
      fprintf('Centering (mean position = %d)\n',meanpos);
      initpos = meanpos - floor((N - 1) / 2);
      finpos = initpos + N - 1;

      finc1 = [];
      for i = initpos:finpos
         if (i < 1 || i > N)
            finc1 = [finc1 meanc1backinde]; 
         else
            finc1 = [finc1 c1(:, i)]; 
         end
      end
      c1 = finc1;
      prevmeanpos = meanpos;
   end

   % No interactions (base independent model)
   hinde=log(c1);
   Savehindepwm{Count}=hinde;

   % Now we find nmers probas to find the sequences predicted by the
   % refined model
   energymodel = energy(nmersint,hinde,Jinde);
   scoresites = energymodel - energybackinde;


   % we sort according to the enrichment compare to background
   [sproba,i]=sort(scoresites,'descend');
   for nmax=1:length(i)
      if length(unique(enhs(i(1:nmax)))) > Nhalf
         Nmax=nmax;
         break;
      end
   end

   seqs=nmers(i(1:Nmax),:);
   Saveseqsinde{Count}=seqs;

   diffsites = size(setdiff(seqs, prevseqs, 'rows'),1);
   fprintf('Iteration %d: %d site(s) differ\n', Count, diffsites);

   Savediffsitespwm(Count)=diffsites;

   Count=Count+1;

   if diffsites==0 && length(seqs) == length(prevseqs)
      fprintf('Converged!\n');
      break
   end
end

fprintf('Saving to %s\n',filePWM);
save(filePWM, 'hinde', 'Jinde', 'seqs');

% uncomment to show logos of the initial and final PWMs (with uniform
% background frequencies):
% seqlogo(initseqs)
% seqlogo(seqs)

fprintf('Exit normally!\n');
