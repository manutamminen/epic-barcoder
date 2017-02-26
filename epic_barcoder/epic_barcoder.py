import os
import subprocess
import string
import random
import time
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


array_dict = {'lsf': '''#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem={}]"
#BSUB -J {}[1-{}]
#BSUB -o wrf.%I
#BSUB -e wrf.%I
#BSUB -W {}

cd {}

name=$(sed -n "$LSB_JOBINDEX"p {})

{}
''', 'slurm': '''#!/bin/bash
#SBATCH --mem-per-cpu={}
#SBATCH -J {}
#SBATCH --array=1-{}
#SBATCH -t {}
#SBATCH -o array_job_out_%j.txt
#SBATCH -e array_job_err_%j.txt
#SBATCH -n 1
#SBATCH -p serial

cd {}

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p {})

{}
'''}


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


def get_seed_dict(uc_file):
    seed_dict = {}
    with open(uc_file) as f:
        for line in f:
            if line.split()[0] == "H":
                seq_id = line.split()[8]
                seed_id = line.split()[9]
                seed_dict[seq_id] = seed_id
            if line.split()[0] == "S":
                seq_id = line.split()[8]
                seed_dict[seq_id] = seq_id
    return seed_dict


def add_otus_to_fasta(seq_file, uc_file, output_file):
    seeds = get_seed_dict(uc_file)
    seq_acc = []
    for seq_id, seq in ep.read_fasta(seq_file):
        short_id = seq_id[1:].split()[0]
        seed_id = seeds[short_id]
        new_seq_id = "{} OTU={}".format(seq_id.strip(), seed_id)
        seq_acc.append([new_seq_id, seq])
    ep.write_fasta(seq_acc, output_file)


def generate_id(size=8):
    """ Generate random sequences of characters for temporary file names.
    """
    chars = string.ascii_uppercase + string.digits
    return ''.join(random.choice(chars) for _ in range(size))


def make_split_seq_dict(seq_iter, no_splits):
    seq_list = list(seq_iter)
    seq_len = len(seq_list)
    chunk_size = int(seq_len/no_splits)
    split_dict = defaultdict(list)
    for ix, (seq_id, seq) in enumerate(seq_list):
        if ix % chunk_size == 0:
            chunk_id = generate_id()
        split_dict[chunk_id].append([seq_id, seq])
    return split_dict


def split_seqs(seq_file, no_splits):
    seqs = ep.read_fasta(seq_file)
    split_dict = make_split_seq_dict(seqs, no_splits)
    for key, val in split_dict.items():
        seq_name = key + "_tmp.fasta"
        ep.write_fasta(val, seq_name)
    return list(split_dict.keys())


def make_array_job(seqs, batch_command, post_command=None, no_splits=1000, scheduler='slurm', memory=2048, run_time='02:00', cleanup=True):
    job_name = generate_id()
    user = subprocess.check_output('whoami', universal_newlines=True).strip()
    namelist = job_name + "_tmp.namelist"
    seq_ids = split_seqs(seqs, no_splits)
    job_no = len(seq_ids)
    home_dir = os.getcwd()
    array = array_dict[scheduler].format(memory, job_name, job_no,
                                         run_time, home_dir, namelist,
                                         batch_command)
    array_file_name = generate_id() + "_tmp.sh"
    with open(array_file_name, "w") as f:
        for line in array:
            f.write(line)
    with open(namelist, "w") as f:
        for item in seq_ids:
            f.write(item + "_tmp.fasta\n")
    if scheduler == 'slurm':
        subprocess.call(['sbatch', array_file_name])
        print("A total of {} jobs.".format(job_no))
        time.sleep(10)
        while True:
            jobs = subprocess.check_output(['squeue', '-u', user],
                                           universal_newlines=True).split("\n")
            if len(jobs) == 2:
                break
            print("{} jobs left".format(len(jobs) - 2))
            time.sleep(5)
    if post_command:
        print("Executing the post-batch command.")
        subprocess.call(post_command.split(" "))
    if cleanup:
        print("Cleaning up.")
        [os.remove(tmp_file) for tmp_file in os.listdir() if "tmp" in tmp_file]
        [os.remove(tmp_file) for tmp_file in os.listdir() if tmp_file.endswith(".txt")]
    print("Done!")


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
