% Defines main input parameters

%% The motif of interest (in this example, Twist in Drosophila)
if (~exist('MAT','var'))
   MAT = 'twi';
   fprintf('Motif: %s\n',MAT);
end

%% !!runMeFirst.sh has to be executed for load_paths file to be created
if (~exist('load_paths','file'))
    fprintf('Pease execute runMeFirst.sh to create load_paths.m file.\n');
    return
end
load_paths;

%% The result directories
outdirPWM = strcat(RESDIR,'/RES-PWM');
filePWM = strcat(outdirPWM,'/',MAT,'.mat');
filePWMback = strcat(outdirPWM,'/',MAT,'-background.mat');

outdirPIM = strcat(RESDIR, '/RES-PIM');
filePIM = strcat(outdirPIM,'/',MAT,'.mat');
filePIMback = strcat(outdirPIM,'/',MAT,'-background.mat');

outdirNNM = strcat(RESDIR, '/RES-NNM');
fileNNM = strcat(outdirNNM,'/',MAT,'.mat');
fileNNMback = strcat(outdirNNM,'/',MAT,'-background.mat');

outdirMIX = strcat(RESDIR, '/RES-MIX');
fileMIX = strcat(outdirMIX,'/',MAT,'.mat');

%% per default, do not restrict to nearest neighbors
if (~exist('nearest_neighbors','var'))
    nearest_neighbors = false;
end

%% DNA alphabet size
if (~exist('q','var'))
    q = 4;
end

%% N-mers size
if (~exist('N','var'))
    N = 12;
    fprintf('N=%d\n',N);
end

%% pseudo-count used in the study
if(~exist('pseudocount','var'));
    pseudocount=1;
    fprintf('pseudocount=%d\n',pseudocount);
end



%% a file giving the N-mers (here N=12) to be searched for binding sites.
% For each ChIPseq peak, a unique identifier was assigned (first column)
% and all N-mers from both strands were considered (second column).
if (~exist('filenmers','var'))
    filenmers = strcat(DATADIR,'/nmers-chips/',MAT,'.dat');
    
    fprintf('Loading ChIPseq N-mers: %s...\n',filenmers);
    fid = fopen(filenmers);
    data = textscan(fid,'%d%s');
    fclose(fid);
    
    nmers = cell2mat(data{2});
    nmersint = double(nt2int(nmers));
    enhs = data{1};
    fprintf('Loaded %d peaks.\n',length(unique(enhs)));
    
    Mchip = size(nmersint,1);
    
    % Nhalf serves as a cutoff to stop learning procedure. It is the number
    % of ChIP peaks that are detected by the model (True Positives), i.e
    % that contain at least one predicted binding site.
    Nhalf = length(unique(enhs))/2; 
end


%% The initial PWM (rows correspond to A,T,C,G frequencies at the different positions of the site)
if (~exist('filemot','var'))
    filemot = strcat(DATADIR,'/motinits/',MAT,'.dat');
    fprintf('Initial PWM file: %s\n',filemot); 
    
    c1init = importdata(filemot);
    c1init = c1init([1,3,4,2],:); % ATCG to ACGT!
    lenc1 = size(c1init, 2);
    % avoid problems with c1=0
    for i=1:lenc1
        c1init(:,i) = (c1init(:,i) + 1e-6) / (sum(c1init(:,i)) + q * 1e-6);
    end
    
end
