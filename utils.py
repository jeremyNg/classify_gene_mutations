from collections import namedtuple

def classify_insertions_deletion(normal, mutant):
    """
    classify insertions/deletions into frameshift or in-frame
    this is done by getting the remainder of difference in length following
    division by 3. If the remainder is zero, it is in-frame
    """
    difference_in_length = len(mutant) - len(normal)
    if difference_in_length % 3:
        # Python treats all non-zero values as True, hence this will be
        # executed when we have a remainder following division by 3
        if len(mutant) > len(normal):
            return "frameshift insertion"
        return "frameshift deletion"
    else:
        if len(mutant) > len(normal):
            return "in-frame insertion"
        return "in-frame deletion"


def map_codon_to_amino_acid(codon):
    """
    converts codon (3 bases) to amino acid
    if the codon is not found in the codon dictionary (the hash map), it will be
    treated as if the codon is a stop codon (*)
    """
    codon_map = {"AAA":"A",
                "TTT":"C",
                "ATG":"M",
                "CTA":"G",
                "CTG":"G",
                "AAT":"W"}
    return  codon_map.get(codon, "*")


def map_amino_acid_to_property(amino_acid):
    """
    function to tell us if a missense substitution results in a replacement of an
    amino acid is chemically similar to the original. We will return None if the
    chemical property cannot be found in our dictionary
    """

    amino_acid_map = {"A":"uncharged",
                       "G": "uncharged",
                       "C": "polar",
                       "M": "polar",
                       "W": "hydrophobic"}
    return amino_acid_map.get(amino_acid, None)

def translate_sequence(sequence):
    """
    translates sequence to amino acid
    """

    amino_acid_sequence = ""

    for start in range(int(len(sequence)/3)):
        try:
            first = start * 3
            end = first + 3
            codon = sequence[first:end]
            amino_acid = map_codon_to_amino_acid(codon)

            if amino_acid == "*":
                """
                we will not add the stop codon, preferring instead to return
                the amino acid sequence as is. this will allow us to immediately
                know when a SNP is a missense one since the amino acid sequence
                will be shorter than the normal one
                """
                return amino_acid_sequence
            else:
                amino_acid_sequence += amino_acid
        except IndexError:
            """
            if the length of provided sequence is not a multiple of 3, we will ignore
            the remaining bases and not attempt to translate it, returning only up to the
            last amino acid added
            """
            return amino_acid_sequence
    return amino_acid_sequence

def align_sequences(normal, mutant):
    """
    aligns two amino acid and identifies locations where they differ
    returns a list of named tuples
    """

    Mutations = namedtuple('mutations', 'normal,mutant')
    mutations_in_sequence = []
    for aa_position, normal_aa in enumerate(normal):
        mutant_aa = mutant[aa_position]
        if mutant_aa != normal_aa:
            mutations_in_sequence.append(Mutations(normal_aa, mutant_aa))
    return mutations_in_sequence


def classify_mutation(normal, mutant):
    """
    main function to classify mutations
    """

    if len(normal)!=len(mutant):
        # If the two lengths differ at nucleotide level, then it will be
        # an insertion/deletion, and we can simply classify as in-frame or
        # frameshift
        return classify_insertions_deletion(normal, mutant)
    else:
        # when the nucleotide lengths are the same, a substitution is implied.
        # hence we need to translate the sequence in order to determine
        # the consequence of the substitution
        normal_sequence = translate_sequence(normal)
        mutant_sequence = translate_sequence(mutant)

        if mutant_sequence == normal_sequence:
            return "Synonymous mutation"
        else:
            if len(mutant_sequence) == len(normal_sequence):
                # Mis-sense substitution
                amino_acid_changes = align_sequences(normal_sequence, mutant_sequence)
                for aa_change in amino_acid_changes:
                    normal_property = map_amino_acid_to_property(aa_change.normal)
                    mutant_property = map_amino_acid_to_property(aa_change.mutant)
                    if normal_property!=mutant_property:
                        return "Missense substitution led to replacement with different properties"
                return "Missense substitution led to replacement with similar properties"
            elif len(mutant_sequence) < len(normal_sequence):
            # Since we terminate translation upon encountering a stop codon,
            # the mutant sequence will be shorter; hence, if there is a gain of
            # stop codon, this block is executed.
                return "Nonsense mutation"
            else:
                # This block is executed only if the amino acid sequence length
                # of the mutant sequence is longer than the normal sequence,
                # implying a loss of stop codon
                return "Loss of stop codon detected following substitution"
