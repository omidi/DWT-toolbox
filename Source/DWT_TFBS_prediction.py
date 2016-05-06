import os
import sys
import subprocess
import math
import re
import textwrap
import numpy as np
from itertools import combinations, product
from string import upper

def arguments():
    import argparse
    parser = argparse.ArgumentParser(description=textwrap.fill(textwrap.dedent("""
    Given a DNA FASTA file and a DWT file, this code predicts TF binding sites.
    """), width=70, break_long_words=True, break_on_hyphens=True),
                                     prog='fitting_DWT_model.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', dest='fasta', action='store',
                        type=str, required=True, metavar='\b',
                        help="""The DNA FASTA file""")
    parser.add_argument('-d', '--dwt', dest='dwt', action='store',
                        type=str, required=True, metavar='\b',
                        help="""The DWT flat file""")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store',
                        type=float, required=False, default=.5, metavar='\b',
                        help="""Optional cutoff over binding posterior (Default 0.5)""")
    parser.add_argument('-o', '--output', dest='output', action='store',
                        type=str, required=True, metavar='\b',
                        help="""Output file name""")
    parser.add_argument('-b', '--with_bg', action='store_true', default=False,
                        help="""Whether the nucleotide background frequency is set according
                        to the DNA sequences. If this option not used, a uniform background
                        frequency will be used.""", dest='with_background')
    args = parser.parse_args()
    if args.cutoff > 1 or args.cutoff < 0:
        print '\nMinimum posterior cutoff is invalid! \n'
        exit()
    return args


def run(cmd):
    proc = subprocess.Popen (cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    if proc.poll() > 0:
        return -1
    return 0



def background_model(infile):
    """
    As input, this method takes a filename which contains the sequences in FASTA format
    and returns the sequences' base frequencies in a form of a dictionary.
    """
    test = re.compile ("^>", re.IGNORECASE)
    reverse = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    freq = {'A':0, 'C':0, 'G':0, 'T':0}
    for line in open(infile):
        if test.search(line):
            continue
        else:
            for b in line.rstrip():
                try:
                    freq[b] = freq[b] + 1
                    freq[reverse[b]] = freq[reverse[b]] + 1
                except KeyError:
                    continue
    normalization = float(sum(freq.values()))
    for b in freq.keys():
        freq[b] = freq[b]/normalization
    return freq



def make_parameter_file(args):
    fname = args.dwt + "_param_file"
    with open(fname, 'w') as outf:
        outf.write('\n'.join([
            "#background",
            "A\t%0.4f" % args.background["A"],
            "C\t%0.4f" % args.background["C"],
            "G\t%0.4f" % args.background["G"],
            "T\t%0.4f" % args.background["T"],
            "",
            "#output_file",
            args.output,
        ]))
    return fname


def pr(s_i, s_0):
    return np.exp(s_i - s_0) / (1 + np.exp(s_i - s_0))


def fit_prior(S):
    def likelihood(S, s_0):
        sum = 0.
        for s in S:
            sum += math.log(1 + math.exp(s - s_0))
        sum -= len(S)*math.log(1+math.exp(-1*s_0))
        return sum
    s0_range = [-10., 30.]
    s_0 = s0_range[0]
    fitted_s_0 = s_0
    max = likelihood(S,s_0)
    for i in xrange(int(math.ceil((s0_range[1]-s0_range[0])/0.1))):
        s_0 += .1
        tmp = likelihood(S,s_0)
        if tmp > max:
            max = tmp
            fitted_s_0 = s_0
    return fitted_s_0


def main():
    args = arguments()
    if args.with_background:
        args.background = background_model(args.fasta)
    else:
        args.background = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    param_file = make_parameter_file(args)

    source_dir = os.path.dirname(sys.argv[0])
    prog = os.path.join(source_dir, "DWT_TFBS_prediction")
    cmd = " ".join([
        prog,
        args.dwt,
        args.fasta,
        param_file
    ])
    if run(cmd) == -1:
        print '\nError in running the following command:'
        print cmd + '\n'
        print 'Program stops! \n'
        exit()

    with open(args.output) as inf:
        output = inf.readlines()

    scores = [float(line.split()[-1]) for line in output]
    prior = fit_prior(scores)

    with open(args.output, 'w') as outf:
        outf.write('\t'.join([
                    "seq",
                    "start",
                    "end",
                    "strand",
                    "TF",
                    "site",
                    "DWT_score",
                    "posterior\n"
                    ]))
        for line in output:
            score = float(line.split()[-1])
            posterior =  pr(score, prior)
            if posterior >= args.cutoff:
                outf.write( '\t'.join([
                    line.rstrip(),
                    '%0.6f\n' % posterior,
                    ]))

    os.system("rm %s" % param_file)

if __name__ == '__main__':
    main()