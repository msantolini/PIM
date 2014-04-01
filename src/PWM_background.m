% Computes the PWM background model given all ChIPseq nmers
load_params;
[~,~] = mkdir(outdirPWM);

% c1back is the frequency matrix of single nucleotide counts
c1back = (histc(nmersint,1:q)+pseudocount)/(size(nmersint,1)+q*pseudocount); 
hindeback = log(c1back);
Jindeback = zeros(q,q,N,N);

meanc1backinde = mean(c1back,2);
energybackinde = energy(nmersint,hindeback,Jindeback);
probabackinde = exp(energybackinde);

fprintf('Saving to %s\n',filePWMback);
save(filePWMback, 'c1back', 'hindeback', 'Jindeback', 'meanc1backinde','energybackinde');