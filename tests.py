from collections import namedtuple
import unittest

from utils import *

class TestUtils(unittest.TestCase):
    # Test cases for insertion and deletions
    def test_frameshift_deletion(self):
        normal = "ATCTTTACG"
        mutant = "ATCTTTCG"
        self.assertEqual(classify_mutation(normal, mutant), "frameshift deletion")

    def test_inframe_deletion(self):
        normal = "ATCTTTACG"
        mutant = "ATCTTT"
        self.assertEqual(classify_mutation(normal, mutant), "in-frame deletion")

    def test_frameshift_insertion(self):
        normal = "ATCTTTACG"
        mutant = "ATCTTTAAAAC"
        self.assertEqual(classify_mutation(normal, mutant), "frameshift insertion")

    def test_inframe_insertion(self):
        normal = "ATCTTTACG"
        mutant = "ATCTTTAAAACG"
        self.assertEqual(classify_mutation(normal, mutant), "in-frame insertion")

    # Test cases for translation sub-routine
    def test_translation(self):
        sequence = "AAATTTATGCTA"
        self.assertEqual(translate_sequence(sequence), "ACMG")

    def test_translation_stop(self):
        sequence = "AAATTTATGCCG" # Since CCG is not in our codon table, it will be treated as a stop
        self.assertEqual(translate_sequence(sequence), "ACM")

    def test_translation_incomplete(self):
        """
        test case for an improper sequence
        (incomplete reading frame should lead to one amino acid less)
        """
        sequence = "AAATTTATGCT"
        self.assertEqual(translate_sequence(sequence), "ACM")

    # test case for sequence alignment and other sub-routines
    def test_align_sequence(self):
        normal = "ATG"
        mutant = "ACG"
        Mutations = namedtuple('mutations', 'normal,mutant')
        result = [Mutations("T", "C")]
        self.assertListEqual(align_sequences(normal, mutant), result)

    def test_map_amino_acid_to_property(self):
        Mutations = namedtuple('mutations', 'normal,mutant')
        self.assertEqual(map_amino_acid_to_property("A"), "uncharged")

    # Test case for SNPs
    def test_nonsense_mutation(self):
        # Since CCG is not in our codon table,
        # it will be treated as a stop and hence should be classified as nonsense
        normal = "AAATTTATGCTA"
        mutant = "AAATTTATGCCG"
        self.assertEqual(classify_mutation(normal, mutant), "Nonsense mutation")

    def test_synonymous_mutation(self):
        # Since CCG is not in our codon table,
        # it will be treated as a stop and hence should be classified as nonsense
        normal = "AAATTTATGCTA"
        mutant = "AAATTTATGCTG"
        self.assertEqual(classify_mutation(normal, mutant), "Synonymous mutation")

    def test_missense_mutation_different_properties(self):
        # Replacement of CTA to TTT leads to change from G (uncharged) to C (polar) at amino acid level
        normal = "AAATTTATGCTA"
        mutant = "AAATTTATGTTT"
        self.assertEqual(classify_mutation(normal, mutant),
                        "Missense substitution led to replacement with different properties")

    def test_missense_mutation_similar_properties(self):
        # Replacement of CTA to TTT leads to change from G to A at last codon
        # Both amino acids are of similar property
        normal = "AAATTTATGCTA"
        mutant = "AAATTTATGAAA"
        self.assertEqual(classify_mutation(normal, mutant),
                        "Missense substitution led to replacement with similar properties")

    def test_stop_lost(self):
        # Test case for stop lost. CCG is used as the stop codon in the normal sequence
        normal = "AAATTTATGCCG"
        mutant = "AAATTTATGCTA"
        self.assertEqual(classify_mutation(normal, mutant),
                        "Loss of stop codon detected following substitution")

if __name__ == '__main__':
    unittest.main()
