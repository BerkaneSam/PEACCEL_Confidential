import random
from urllib import request
import os
import argparse

aat = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


def get_arguments():
    parser = argparse.ArgumentParser(description="create files to be used in InnovSar")
    parser.add_argument('data', metavar='data', type=str, help="path to pretreated data to be used, sequences must "
                                                               "always be in the last cells")
    parser.add_argument('-m', '--mutpos', type=int, required=True, help="position of mutation in csv file (first cell "
                                                                        "is 0)")
    parser.add_argument('-f', '--fitpos', type=int, required=True, help="position of fitness in csv file (first cell "
                                                                        "is 0")
    parser.add_argument('-ff', '--fasta', type=str, help="Retrieve WT sequence from fastafile, default retrieve WT "
                                                         "sequence from uniprot using data filename, exemple:"
                                                         "P06657_data.csv, in the default case _ is necessary after the"
                                                         "uniprot code in first position")
    parser.add_argument('-wt', '--wildtype', type=str, help="Make use of WT sequence")
    return parser.parse_args()


def check_aa(seq, aas):
    """
    checks if all amino acids are found in a list of known amino acids returns 1 if there is an unknown amino acid
    :param seq: sequence of amino acids
    :param aas: list of known amino acids
    :return: a boolean value
    """
    for aa in seq:
        if aa not in aas:
            return 1
    return 0


def data_retrieve(path, aas, cellmut, cellfit):
    """
    Retrieves data of interest (variants, sequence, fitness, mutation) from dataset and check if the sequence is clean
    and if all those data of interest are available in the dataset
    :param path: path to the dataset
    :param aas: list of amino acids
    :param cellmut: column in which to find the mutation in the dataset
    :param cellfit: column in which to find the fitness in the dataset
    :return: a list of variants each containing a list of the data of interest
    """
    print("retrieving data")
    dataset = []
    with open(path, 'r')as filin:
        for line in filin:
            if not line.startswith('id'):
                sline = line.rstrip().split(',')
                temp_dat = [sline[-1], sline[cellfit], sline[cellmut]]
                if temp_dat[1] != '':
                    if check_aa(temp_dat[0], aas) == 0:
                        dataset.append(temp_dat)
    return dataset


def set_making(dataset):
    """
    Making training set and testing set by shuffling data and separation by a 0.8/0.2 ratio
    :param dataset: list of list of data
    :return: two lists, training set and testing set
    """
    random.shuffle(dataset)
    datacut = int(len(dataset) * 0.8)
    testset = dataset[datacut:]
    trainset = dataset[:datacut]
    print(f"trainset size : {len(trainset)}")
    print(f"testset size : {len(testset)}")
    return trainset, testset


def data_sep(data):
    """
    Retrieving the fitness and sequences into two different lists
    :param data: data to be used
    :return: a list of sequences and a list of fitness
    """
    seq = []
    fit = []
    for intel in data:
        seq.append(intel[0])
        fit.append(intel[1])
    return seq, fit


def sort_key(data):
    """
    used to select a sorting key
    :param data: data to be sorted
    :return: key
    """
    return float(data[1])


def mut_comb_s(dataset):
    """
    Retrieving the mutations with the best fitness to make a file by sorting the dataset variants by fitness
    :param dataset: list of list of data
    :return: dict of 10 best mutations
    """
    sorted_fit = sorted(dataset, key=sort_key, reverse=True)
    muta = []
    for data in sorted_fit:
        muta.append(data[2])
        if len(muta) >= 10:
            return muta
    return muta


def mut_writing_s(mut_dict, direc):
    """
    Making mutation file
    :param mut_dict: list of mutations
    :param direc: part of filename
    :return: makes a file
    """
    with open(f"{direc}/mut_comb.txt", 'w')as filout:
        for mut in mut_dict:
            filout.write(f"{mut}\n")


def seq_retrieve(unifile):
    """
    Retrieving sequence from a html uniprot webpage
    :param unifile: path to html file
    :return: sequence of amino acids
    """
    seq = ""
    with open(unifile, 'r')as filin:
        for i in range(24):
            next(filin)
        for line in filin:
            if not line.startswith('<'):
                seq = seq + line.rstrip()
            if line.startswith('<'):
                break
    return seq


def retrieve_wt(data_path):
    """
    Retrieve a html uniprot webpage by using the uniprot code of a protein and then retrieves a sequence of amino acids
    in the page, deletes the html page when use of the function is over
    :param data_path: name of a file with formated name
    :return: a sequence of amino acids
    """
    name = data_path.split('_')
    request.urlretrieve(f'https://www.uniprot.org/uniprot/{name[0]}', 'unifile.txt')
    wt_seq = seq_retrieve('unifile.txt')
    os.remove('unifile.txt')
    return wt_seq


def writing_wt(seq, direc):
    """
    Making a file for the sequence of reference
    :param seq: sequence of amino acids
    :param direc: part of filename
    :return: makes a file
    """
    with open(f"{direc}/wt_seq.txt", 'w')as filout:
        filout.write(seq)


def writing_files(data, filename):
    """
    Makes files of training set, testing set, fitness...
    :param data: list to be used
    :param filename: name of file to be made
    :return: makes a file
    """
    with open(filename, 'w')as filout:
        for value in data:
            filout.write(value + "\n")


def making_directory(path):
    """
    Retrieving path of directory where the data are kept and where outputs are to be saved
    :param path: path of dataset
    :return: path to directory
    """
    if '/' in path:
        splpath = path.split('/')
        name = splpath[-1].split('.')
    else:
        name = path.split('.')
    os.mkdir(name[0])
    return name[0]


def retrieve_seq(fastapath):
    """
    retrieves a sequence of amino acids from a fasta file
    :param fastapath: path to the fasta file
    :return: a sequence of amino acids
    """
    sequence = ""
    with open(fastapath, 'r')as filin:
        for line in filin:
            if not line.startswith('>'):
                sequence += line.rstrip()
    return sequence


if __name__ == '__main__':
    print("program launched")
    args = get_arguments()
    print("starting data retrieving...")
    main_data = data_retrieve(args.data, aat, args.mutpos, args.fitpos)
    if args.fasta:
        seq_wt = retrieve_seq(args.fasta)
    else:
        seq_wt = retrieve_wt(args.data)
    print("retrieving done")
    print("making training set and testing set...")
    train_set, test_set = set_making(main_data)
    train_seq, train_fit = data_sep(train_set)
    test_seq, test_fit = data_sep(test_set)
    print("set making done")
    print("retrieving mutations combinations...")
    mutac = mut_comb_s(main_data)
    print("mutations retrieved")
    print("writing files...")
    direct = making_directory(args.data)
    writing_files(train_seq, f"{direct}/trainset_seq.txt")
    writing_files(train_fit, f"{direct}/trainset_fit.txt")
    writing_files(test_seq, f"{direct}/testset_seq.txt")
    writing_files(test_fit, f"{direct}/testset_fit.txt")
    writing_wt(seq_wt, direct)
    mut_writing_s(mutac, direct)
    print("writing done")
    print("program done")
