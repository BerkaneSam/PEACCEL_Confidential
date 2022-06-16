import argparse


aat = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


def get_arguments():
    parser = argparse.ArgumentParser(description="Checks if a fasta sequence corresponds to a dataset with mutation "
                                                 "informations found in separate columns")
    parser.add_argument('data', metavar='data', type=str, help="the dataset to be used")
    parser.add_argument('-f', '--fasta', required=True, help="fasta file to be used")
    parser.add_argument('-aa', '--originaa', required=True, help="position of original aa in file")
    parser.add_argument('-pos', '--position', required=True, help="position of aa's position in sequence in the file")
    parser.add_argument('-j', '--jumpl', type=int, help="number of lines to jump in file")
    return parser.parse_args()


def retrieve_seq(fastapath):
    """
    Retrieves sequence of amino acids from fasta file
    :param fastapath: path to fasta file
    :return: a sequence of amino acids
    """
    sequence = ""
    with open(fastapath, 'r')as filin:
        for line in filin:
            if not line.startswith('>'):
                sequence += line.rstrip()
    return sequence


def retrieve_pos(datapath, aa, trueaapos, posp, jump=0):
    """
    Retrieves original aa and position in the sequence from dataset(from mutation) for each variants and adds them to a
    list and checks if only natural amino acids are present
    :param datapath: path to dataset file
    :param aa: list of natural amino acids
    :param trueaapos: column where to find original amino acid for each variants
    :param posp: column where to find position of amino acids in the sequence for each variants
    :param jump: number of line to jump to reach variants
    :return: returns a list of original amino acids and their position in the sequence (aapos)
    """
    true_pos = []
    with open(datapath, 'r')as filin:
        for i in range(jump):
            next(filin)
        for line in filin:
            sline = line.rstrip().split(',')
            mutpos = sline[int(trueaapos)] + sline[int(posp)]
            if mutpos[0] in aa:
                true_pos.append(mutpos)
    return true_pos


def check_correspondence(fastaseq, aapos, namefile):
    """
    Compares the original amino acid found in mutation data from dataset to the reference sequence retrieved to check
    if they corresponds and write a log file to see the position where there are differences if there is any
    :param fastaseq: reference sequence
    :param aapos: list of original amino acids and their position in the sequence
    :param namefile: new file name
    :return: a boolean value as well as writing a log file
    """
    check = 0
    with open(f'{namefile}_checking_log.txt', 'w')as filout:
        for posaa in aapos:
            try:
                pos = int(posaa[1:]) - 1
            except ValueError:
                continue
            trueaa = posaa[0]
            try:
                if fastaseq[pos] == trueaa:
                    filout.write(f'{posaa} fasta aa = {fastaseq[pos]}  true aa = {trueaa}  checked -> True\n')
                else:
                    check = 1
                    filout.write(f"{posaa} fasta aa = {fastaseq[pos]}  true aa = {trueaa}  checked -> False -----------"
                                 f"-----------\n")
            except IndexError:
                filout.write(f"{posaa} Not found in sequence's range|Sequence too short(most likely) checked -> False"
                             f"------------------\n")
                check = 1
        if check == 0:
            filout.write("\nFull sequence checked\nResult : Sequence corresponds to dataset")
        else:
            filout.write("\nFull sequence checked\nResult : Sequence does not correspond to dataset")
    return check


def get_name(name):
    """
    Treats file name to retrieve the name without directory and file format
    :param name: path to file
    :return: striped file name
    """
    if "/" in name:
        partname = name.split("/")
        truename = partname[-1].split(".")
    else:
        truename = name.split(".")
    return truename[0]


def get_logname(dataname, fastaname):
    """
    Makes a new filename for logs based on dataset name and fasta name
    :param dataname: path to dataset
    :param fastaname: path to fasta file
    :return: new log file name
    """
    data_truename = get_name(dataname)
    fasta_truename = get_name(fastaname)
    realname = f"{data_truename}_{fasta_truename}"
    return realname


if __name__ == '__main__':
    print("program launched")
    args = get_arguments()
    log_name = get_logname(args.data, args.fasta)
    print("arguments retrieved")
    print("retrieving positions from dataset...")
    if args.jumpl:
        data = retrieve_pos(args.data, aat, args.originaa, args.position, args.jumpl)
    else:
        data = retrieve_pos(args.data, aat, args.originaa, args.position)
    print("position retrieved")
    print("retrieving sequence...")
    seq = retrieve_seq(args.fasta)
    print("sequence retrieved")
    print("checking sequence...")
    result = check_correspondence(seq, data, log_name)
    print("sequence checked and log file written")
    if result == 0:
        print("Sequence validated")
    else:
        print("Sequence not validated")
    print("program done")
