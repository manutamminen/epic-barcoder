import os
from itertools import zip_longest
from collections import defaultdict, Counter
import epride as ep
import pandas as pd


bridges_dict = {"16S": "GWATTACCGCGGCKGCTGCATCTTCTCCAAATGGGTCATGATC",
                "18S": "AAGAACGGCCATGCACCACCACATCTTCTCCAAATGGGTCATGATC",
                "narG2": "ACCGACATGCCGWTSCTGGTCATCTTCTCCAAATGGGTCATGATC",
                "nosZ2": "AACAAGTTCTCSAAGGACCGCATCTTCTCCAAATGGGTCATGATC",
                "nosZ3": "CTCMAAGGACCGGTTCMTSMCATCTTCTCCAAATGGGTCATGATC",
                "norB2": "GCACCGGYCACCAYTAYTWCCATCTTCTCCAAATGGGTCATGATC"}

size_filter_dict = {'narG2': 100, 'norB2': 90, 'nosZ3': 100, '18S': 70, '16S': 90, 'nosZ2': 100}


def move_barcodes_and_type_to_fasta_id(bc_seq, bridge_dict):
    bridge_dict = {key: ep.expand_primers(ep.reverse_complement(val))
                   for key, val in bridge_dict.items()}
    for seq_id, seq in bc_seq:
        for bridge_id, bridges in bridge_dict.items():
            for bridge in bridges:
                if bridge in seq:
                    bc, rest = seq.split(bridge)
                    if len(bc) == 20:
                        seq_id = "{} barcode={} sequence_type={}".format(seq_id.strip(),
                                                                        bc, bridge_id)
                        yield([seq_id, rest])


def process_barcode_info(input_seq_file, output_seq_file, bridge_dict):
    seqs = ep.read_fasta(input_seq_file)
    fasta_iter = move_barcodes_and_type_to_fasta_id(seqs, bridge_dict)
    ep.write_fasta(fasta_iter, output_seq_file)


def get_len_distr(seqs):
    len_distr_dict = defaultdict(list)
    for seq_id, seq in seqs:
        seq_type_id = seq_id.split("_")[-1]
        len_distr_dict[seq_type_id].append(int(len(seq) / 10) * 10)
    return len_distr_dict


class BCSeq(object):
    def __init__(self, seq, bridges, length_dict=None):
        self.bridges = bridges
        rc_bridges = {key: ep.expand_primers(ep.reverse_complement(val))
                      for key, val in self.bridges.items()}
        self.seq_list = list(ep.read_fasta(seq))
        self.processed_seqs = list(move_barcodes_and_type_to_fasta_id(self.seq_list, rc_bridges))
        self.seq_dict = defaultdict(list)
        for seq_id, seq in self.processed_seqs:
            seq_short_id = seq_id.split("_")[0].split(">")[1]
            self.seq_dict[seq_short_id].append([seq_id, seq])
        self.bc_data = {}
        self.li = LenDist({key: Counter(map(lambda x: int(x / 10) * 10, val)) for key, val in
                           get_len_distr(self.processed_seqs).items()}, length_dict)
        self.length_dict = length_dict

    def __getitem__(self, ix):
        if ix in self.bc_data:
            return self.bc_data[ix]
        else:
            self.bc_data[ix] = BCData(self.seq_dict[ix], self.length_dict)
            return self.bc_data[ix]

    def __repr__(self):
        repr_lst = ["Sequence group: {:15} number of sequences: {}".format(seq_name, len(seqs))
                    for seq_name, seqs in self.seq_dict.items()]
        repr_str = "\n".join(repr_lst)
        return repr_str


class BCData(object):
    def __init__(self, seq, length_dict=None):
        self.seq_list = seq
        self.li = LenDist({key: Counter(map(lambda x: int(x / 10) * 10, val))
                           for key, val in get_len_distr(self.seq_list).items()}, length_dict)
        self.seq_filt = []
        self.length_dict = length_dict
        if self.length_dict:
            self.filter_seqs()
            bc_type_table = pd.DataFrame([seq_id.strip().split("_")[2:]
                                          for seq_id, seq in self.seq_filt],
                                         columns=["Barcode", "Seq_type"])
        else:
            bc_type_table = pd.DataFrame([seq_id.strip().split("_")[2:]
                                          for seq_id, seq in self.seq_list],
                                         columns=["Barcode", "Seq_type"])
        self.type_table = bc_type_table.groupby(["Barcode", "Seq_type"])\
                                       .size().reset_index().rename(columns={0: "Count"})
        self.bc_table = pd.pivot_table(self.type_table, values="Count", index="Barcode",
                                       columns="Seq_type", fill_value=0)

    def __repr__(self):
        abund_table = self.bc_table.groupby(list(self.bc_table.columns))\
            .size().reset_index().sort_values(0, ascending=False)\
                                 .rename(columns={0: 'Count'}).head(20)
        abund_table.index = range(len(abund_table.index))
        repr_string = abund_table.to_string()
        return repr_string

    def filter_seqs(self):
        for seq_id, seq in self.seq_list:
            seq_type = seq_id.split("_")[-1]
            min_len = self.length_dict[seq_type]
            if len(seq) >= min_len:
                self.seq_filt.append([seq_id, seq])


class LenDist(object):
    def __init__(self, distr_info, length_dict=None):
        self.distr_info = distr_info
        self.length_dict = length_dict

    def process_repr(self):
        for item, val in self.distr_info.items():
            yield([item, ''])
            yield(['LI', 'Count'])
            count_per_item = sorted([[len_interval, count]
                                     for len_interval, count in val.items()],
                                    key=lambda x: x[0])
            for count_item, count in count_per_item:
                try:
                    if self.length_dict[item] == count_item:
                        yield(["---", "---"])
                    yield([count_item, count])
                except TypeError:
                    yield([count_item, count])


    def koe_repr(self):
        tmp = self.process_repr()
        return(list(tmp))
        # acc = []
        # for i in zip(*tmp):
        #     acc.append(i)
        # return acc
        # for i in zip(*self.process_repr()):
        #     acc.append("\t\t".join(i))
        # repr_str = "\n".join(acc)
        # return repr_str

    # def __repr__(self):
    #     # return str(self.process_repr())
    #     acc = []
    #     for i in zip(*self.process_repr()):
    #         acc.append("\t\t".join(i))
    #     repr_str = "\n".join(acc)
    #     return repr_str



# os.chdir("/Users/tamminma/gits/epic_barcoder/epic_barcoder/")

# test_fasta = list(ep.read_fasta("../test_seqs.fasta"))

# a = BCSeq("../test_seqs.fasta", bridges_dict, size_filter_dict)

# print(a)

# print(a['WW'])

# print(a.li.koe_repr())

# b = a.li.koe_repr()
