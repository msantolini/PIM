function seqsbin = nt2bin(seqs)

seqsbin=[];
for i=1:size(seqs,1)
    seqbin=[];
    for j=1:size(seqs,2)
        if (seqs(i,j) == 'A') 
            seqbin=[seqbin 1 0 0 0];
        elseif (seqs(i,j) == 'C') 
            seqbin=[seqbin 0 1 0 0];
        elseif (seqs(i,j) == 'G') 
            seqbin=[seqbin 0 0 1 0];
        elseif (seqs(i,j) == 'T') 
            seqbin=[seqbin 0 0 0 1];
        else
            seqbin=[seqbin 0 0 0 0];
        end
    end
    seqsbin=[seqsbin; seqbin];
end