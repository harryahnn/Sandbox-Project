


function aminoAcidSequence = dNAtoAminoAcid(dnaSequence)
% Step 1: DNA to mRNA Conversion (T -> A, A -> U, C -> G, G -> C)
mRNAsequence = strrep(dnasequence, 'T', 'P'); % Replace T with placeholder so it doesn't overwrite in the next line
mRNAsequence = strrep(mRNAsequence, 'A', 'U'); % Replace A with U
mRNAsequence = strrep(mRNAsequence, 'P', 'A'); % Replaces placeholder with A, which was originally T
mRNAsequence = strrep(mRNAsequence, 'C', 'P'); % Replace C with placeholder
mRNAsequence = strrep(mRNAsequence, 'G', 'C'); % Replace G with C
mRNAsequence = strrep(mRNAsequence, 'P', 'G'); % Replace placeholder with G

% Step 2: Split mRNA into codons (sets of 3 nucleotides)
codonChain = {};
for i = 1:3:length(mRNAsequence)-2
    codonChain{end+1} = mRNAsequence(i:i+2); % Create codons
end

% Step 3: Translate mRNA codons to amino acids using a genetic code table
% Codon to Amino Acid mapping (standard genetic code)
codonTranslated = containers.Map(...
    {'UUU', 'UUC', 'UUA', 'UUG', 'UCU', 'UCC', 'UCA', 'UCG', 'UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG', 'CUU', 'CUC', 'CUA', 'CUG', 'CCU', 'CCC', 'CCA', 'CCG', 'CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG', 'AUU', 'AUC', 'AUA', 'AUG', 'ACU', 'ACC', 'ACA', 'ACG', 'AAU', 'AAC', 'AAA', 'AAG', 'AGU', 'AGC', 'AGA', 'AGG', 'GUU', 'GUC', 'GUA', 'GUG', 'GCU', 'GCC', 'GCA', 'GCG', 'GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG'}, ...
    {'Phe', 'Phe', 'Leu', 'Leu', 'Ser', 'Ser', 'Ser', 'Ser', 'Tyr', 'Tyr', 'Stop', 'Stop', 'Cys', 'Cys', 'Stop', 'Trp', 'Leu', 'Leu', 'Leu', 'Leu', 'Pro', 'Pro', 'Pro', 'Pro', 'His', 'His', 'Gln', 'Gln', 'Arg', 'Arg', 'Arg', 'Arg', 'Ile', 'Ile', 'Ile', 'Met', 'Thr', 'Thr', 'Thr', 'Thr', 'Asn', 'Asn', 'Lys', 'Lys', 'Ser', 'Ser', 'Arg', 'Arg', 'Val', 'Val', 'Val', 'Val', 'Ala', 'Ala', 'Ala', 'Ala', 'Asp', 'Asp', 'Glu', 'Glu', 'Gly', 'Gly', 'Gly', 'Gly'});

% Initialize the amino acid sequence
aminoAcidSequence = '';

% Translate each codon into an amino acid
for i = 1:length(codonChain)
    codon = codonChain{i};
    if isKey(codonTranslated, codon)
        if (i == length(codonChain))
            aminoAcidSequence = [aminoAcidSequence, codonTranslated(codon)];
        else
            aminoAcidSequence = [aminoAcidSequence, codonTranslated(codon), '-'];
        end
    else
        aminoAcidSequence = [aminoAcidSequence, '-', 'X']; % For unknown codons
    end
end
end

function potentialMutation = findMutation(sequence)
    updatedSequence = strsplit(sequence, '-');
    mutationIndices = find(updatedSequence == 'X');
    potentialMutation = mutationIndices;
    if ~isempty(mutationIndices)
        disp('Potential mutations (X) found at the following amino acid positions: ');
        disp(mutationIndices);
    else
        disp('No potential mutations (X) found in this amino acid sequence.');
    end
end
