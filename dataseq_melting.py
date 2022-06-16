import argparse


def get_arguments():
    parser = argparse.ArgumentParser(description="Adds mutated sequence to csv file")
    parser.add_argument('data', metavar='data', type=str, help="the dataset to be used")
    parser.add_argument('-f', '--fasta', required=True, help="fasta file to be used")
    parser.add_argument('-m', '--mutcell', required=True, type=int, help="cell to get mutation in csv file (first cell"
                                                                         "being 0!)")
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


def seq_muting(line, mutpos, sequence):
    """
    Mutating reference sequence based on the mutation of a given variant and checking at the same time if there are no
    mistakes then adding the mutated sequence at the last position of the given line
    :param line: variant line in the dataset containing all data pertaining to the variant
    :param mutpos: column of the mutations
    :param sequence: reference sequence
    :return: new line with the mutated sequence in the last column (added)
    """
    sline = line.rstrip().split(',')
    if ":" not in sline[mutpos]:
        tempseq = list(sequence)
        baseaa = sline[mutpos][0]
        try:
            posaa = int(sline[mutpos][1:-1]) - 1
        except ValueError:
            return line.rstrip() + f",{sequence}\n"
        mutaa = sline[mutpos][-1]
        if tempseq[posaa] == baseaa:
            tempseq[posaa] = mutaa
            templine = line.rstrip() + f",{''.join(tempseq)}\n"
            return templine
    elif ":" in sline[mutpos]:
        tempseq = list(sequence)
        splmut = sline[mutpos].split(':')
        for mut in splmut:
            baseaa = mut[0]
            posaa = int(mut[1:-1]) - 1
            mutaa = mut[-1]
            if tempseq[posaa] == baseaa:
                tempseq[posaa] = mutaa
        return line.rstrip() + f",{''.join(tempseq)}\n"


def filename(data):
    """
    Making new file name
    :param data: name of base file
    :return: new file name
    """
    bdata = data.split('.')
    new_name = bdata[0]+'_seq.'+bdata[1]
    return new_name


def seq_adding(dataset, mutpos, sequence, namef):
    """
    Making new file containing the mutated sequence in the last column by using previous dataset not having sequences
    :param dataset: base dataset
    :param mutpos: column of the mutations
    :param sequence: reference sequence
    :param namef: name of new file
    :return: makes a file
    """
    with open(namef, 'w')as filout:
        with open(dataset, 'r')as filin:
            firstline = next(filin)
            filout.write(firstline.rstrip() + f",sequence\n")
            for line in filin:
                newline = seq_muting(line, mutpos, sequence)
                filout.write(newline)
    print("file wrote")


if __name__ == '__main__':
    print("program launched")
    args = get_arguments()
    nname = filename(args.data)
    print("retrieving WT sequence...")
    wtseq = retrieve_seq(args.fasta)
    print("sequence retrieved")
    print("starting sequence mutations...")
    seq_adding(args.data, args.mutcell, wtseq, nname)
    print("program done")
