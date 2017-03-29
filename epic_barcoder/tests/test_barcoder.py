import unittest
import os
import epride as ep
from .. import epic_barcoder

os.chdir("epic_barcoder/tests/")

class TestUtilities(unittest.TestCase):

    def test_filter_bridge(self):
        fna = list(ep.read_fasta("test.fna"))
        filt = list(epic_barcoder.filter_bridge(fna,
                                                '16S',
                                                'GWATTACCGCGGCKGCTGCATCTTCTCCAAATGGGTCATGATC'))
        self.assertEqual(len(filt), 25)


    def test_filter_reverse(self):
        fna = list(ep.read_fasta("test.fna"))
        filt = list(epic_barcoder.filter_reverse(fna,
                                                 'ATTAGAWACCCBDGTAGTCC'))
        self.assertEqual(len(filt), 34)


    def test_filter_bridge_and_reverse(self):
        fna = list(ep.read_fasta("test.fna"))
        bridge_filt = epic_barcoder.filter_bridge(fna, '16S',
                                                  'GWATTACCGCGGCKGCTGCATCTTCTCCAAATGGGTCATGATC')
        bridge_reverse_filt = epic_barcoder.filter_reverse(bridge_filt,
                                                           'ATTAGAWACCCBDGTAGTCC')
        self.assertEqual(len(list(bridge_reverse_filt)), 13)


    def test_process_barcode_info(self):
        epic_barcoder.process_barcode_info("test.fna", "tmp.fasta", '16S',
                                           'GWATTACCGCGGCKGCTGCATCTTCTCCAAATGGGTCATGATC',
                                           'ATTAGAWACCCBDGTAGTCC')
        tmp_fna = list(ep.read_fasta("tmp.fasta"))
        os.remove("tmp.fasta")
        self.assertEqual([len(seq) for _, seq in tmp_fna],
                         [249, 249, 250, 250, 249, 249, 249, 250, 250, 250, 250, 249, 250])
