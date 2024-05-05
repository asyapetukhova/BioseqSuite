import unittest
import tempfile
import numpy as np
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from Bio import SeqIO

from BioseqSuite import AminoAcidSequence, DNASequence, RNASequence
from custom_random_forest import RandomForestClassifierCustom


class TestRandomForestClassifierCustom(unittest.TestCase):
    def setUp(self):
        self.X, self.y = make_classification(n_samples=1000, n_features=10, n_informative=5, random_state=42)
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size=0.2, random_state=42)
        self.clf = RandomForestClassifierCustom(n_estimators=10, max_depth=5, max_features=3, random_state=42)

    def test_fit(self):
        self.clf.fit(self.X_train, self.y_train)
        self.assertEqual(len(self.clf.trees), self.clf.n_estimators)
        self.assertEqual(len(self.clf.feat_ids_by_tree), self.clf.n_estimators)
        self.assertEqual(len(self.clf.random_states_by_tree), self.clf.n_estimators)

    def test_predict(self):
        self.clf.fit(self.X_train, self.y_train)
        y_pred = self.clf.predict(self.X_test)
        self.assertIsInstance(y_pred, np.ndarray)
        self.assertEqual(y_pred.shape, (self.X_test.shape[0],))

    def test_accuracy(self):
        self.clf.fit(self.X_train, self.y_train)
        y_pred = self.clf.predict(self.X_test)
        accuracy = accuracy_score(self.y_test, y_pred)
        self.assertGreater(accuracy, 0.8)

    def test_invalid_n_estimators(self):
        with self.assertRaises(ValueError):
            RandomForestClassifierCustom(n_estimators=-1, max_depth=5, max_features=3, random_state=42)

class TestSequences(unittest.TestCase):
    def test_amino_acid_sequence(self):
        seq = AminoAcidSequence("ARNDCQEGHILKMFPSTWYV")
        self.assertTrue(seq.is_valid(seq._sequence))

    def test_dna_sequence_complement(self):
        seq = DNASequence("ATCG")
        self.assertEqual(str(seq.complement()), "TAGC")

    def test_dna_sequence_transcribe(self):
        seq = DNASequence("ATCGAT")
        rna_seq = seq.transcribe()
        self.assertIsInstance(rna_seq, RNASequence)
        self.assertEqual(str(rna_seq), "AUCGAU")

class TestFileIO(unittest.TestCase):
    def test_read_write_fasta(self):
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
            tmp.write(">seq1\nARNDCQEGHILKMFPSTWYV\n>seq2\nATCGATCG\n")
            tmp.seek(0)
            sequences = [str(rec.seq) for rec in SeqIO.parse(tmp.name, "fasta")]
            self.assertEqual(sequences, ["ARNDCQEGHILKMFPSTWYV", "ATCGATCG"])

if __name__ == '__main__':
    unittest.main()
