
def arguments():
    import argparse
    parser = argparse.ArgumentParser(description="""Takes a list of binding sites,
    and generates a DWT object.""")

    parser.add_argument('-i', dest='input_file', action='store', required=True, \
                        help="""A list of same size binding sites.
                        Note that the sequences can be weighted by the second column,
                        otherwise each of them will be counted equally (weighted 1.0).""")
    parser.add_argument('-t', dest='TF_name', action='store', required=False, \
                        help='Name of the TF (not mandatory).')
    args = parser.parse_args()
    return args

def pairwise_frequencies(sequences, length):
    convert = {'A':0, 'C':1, 'G':2, 'T':3}    
    dinucleotide_matrix = np.zeros(length*length*16).reshape(length, length, 16)
    for s in sequences:
        for p in combinations(range(length), 2):
            dinucleotide_matrix[p[0]][p[1]][ convert[s[0][p[0]].upper()]*4 + convert[s[0][p[1]].upper()] ] += s[1]
    return dinucleotide_matrix
            

if __name__ == '__main__':
    import numpy as np
    from itertools import combinations, product    
    args = arguments()
    # test if the input file contains two columns or only one
    sequences = []
    f = open(args.input_file)
    if len(f.readline().split('\t')) == 2:
        f.seek(0)
        sequences = [(s.split()[0], float(s.split()[1])) for s in f]
    else:
        f.seek(0)
        sequences = [(s.split()[0], 1.0) for s in f]

    length = len(sequences[0][0])
    dinucleotide_matrix = pairwise_frequencies(sequences, length)
    header = [''.join(i) for i in product('ACGT',repeat=2)]

    if args.TF_name is None:
        print 'NA\tTF'
    else:
        print 'NA\t%s' % (args.TF_name)


    print 'PO1\tPO2\t%s' % ('\t'.join(header))
    for p in combinations(range(length), 2):
        print '%d\t%d\t%s' % (p[0]+1, p[1]+1, '\t'.join(map(str, np.round(dinucleotide_matrix[p[0]][p[1]][:],3) ) ) )

    
    
