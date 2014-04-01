% Computes the NNM background model given all ChIPseq nmers. 
load_params;
[~,~] = mkdir(outdirNNM);

nearest_neighbors = true;
PIM_background;