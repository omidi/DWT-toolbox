
import os
import sys
import subprocess
import math
import re
import numpy as np
from itertools import combinations, product
from string import upper


def arguments():
    import argparse
    parser = argparse.ArgumentParser(description="""
    By providing a DNA FASTA file and a PSWM file (normal weight matrix),
    this codes fit a DWT model, iteratively. The output ia a DWT flat file
    which will be copied in the current working directory or the destination
    directory.
    """, prog='fitting_DWT_model.py', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-w', '--wm', dest='WM',
                      help="""A normal weight matrix file, for formatting consult with the
                      manual.txt file""", \
                      action='store', type=str,  metavar='\b', required=True)
    parser.add_argument('-f' , '--fasta', dest='fasta_file' ,
                      help="""DNA FASTA file""", \
                      action='store', type=str,  metavar='\b', required=True)
    parser.add_argument('-o' , '--outdir', dest='output_dir' ,
                      help="""Optional input for the output directory. if not provided,
                      the current working directory will be used as the output directory.""", \
                      action='store', type=str, required=False,  metavar='\b')
    parser.add_argument('-b', '--with_bg', action='store_true', default=False,
                      help="""Whether the nucleotide background frequency is set according
                      to the provided DNA sequences. But if not used, a uniform background
                      frequency will be used.
                      """, dest='with_background')
    parser.add_argument('-p', '--min_post', action='store', required=False, metavar='\b',
                        help="""Optional minimum posterior cutoff, for selecting TFBS at each
                             round of iteration. (Default 0.5)""",
                        dest="min_post", type=float, default=0.5)
    parser.add_argument('-t', '--tf', action='store', metavar='\b', required=False,
                        help="""Optional TF name, if not given the TF name is assumed to be
                        identical to the name of the PSWM file.
                        """, dest="TF", type=str)
    args = parser.parse_args()
    if args.min_post > 1 or args.min_post < 0:
        print '\nMinimum posterior is invalid! \n'
        exit()
    return args


def mkdir(name):
    proc = subprocess.Popen ('mkdir %s' % name,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    shell=True
                             )
    stdout_value, stderr_value = proc.communicate()
    if proc.poll() > 0:
        if stderr_value.rfind('File exists') > 0:
            # sys.stderr.write ( "\nWarning in creating directory\n")
            # print '\tstderr:', repr(stderr_value.rstrip())
            # print
            # print """The program continues, but it may remove/change file(s) that are
            # already in %s directory""" % name
            None
        else:
             sys.stderr.write ( "\nError in creating directory REGIONS\n" )
             print '\tstderr:', repr(stderr_value.rstrip())
             print 'Program halts!\n'
             sys.exit(-1)
    return True


def run(cmd):
    proc = subprocess.Popen (cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    if proc.poll() > 0:
        return stderr_value.rstrip()
    return None


def runMotevo(inputSequences, paramFile, WM, DirName, program_dir = ""):
    curr_dir = os.getcwd()
    prog = os.path.join(program_dir, 'libs/motevo_1.03/bin/motevo')
    cmd = ' '.join([prog,
                  inputSequences,
                  os.path.join(DirName, paramFile),
                  WM])
    if run(cmd):
        sys.stderr.write ( "\nError in running MotEvo\n")
        print 'Program halts! '
        exit()
    run('rm wms.updated')
    return 0


def motevoParamFile(fileName, motifLength, TFname, args):
    paramFile = open(fileName, 'w')
    if not args.with_background:
        paramFile.write('\n'.join(['Mode WMREF',
                                   'refspecies hg19',
                                   'TREE (hg19:1);',
                                   'wmdiff 0.05',
                                   'EMprior 1',
                                   'priordiff 0.01',
                                   'markovorderBG 0',
                                   'bgprior 0.999',
                                   'bg A 0.25', 'bg T 0.25', 'bg C 0.25', 'bg G 0.25',
                                   'restrictparses 0',
                                   'sitefile %s' % os.path.join(args.output_dir, 'sites'),
                                   'priorfile %s' % os.path.join(args.output_dir, 'priors'),
                                   'minposterior %f' % args.min_post,
                                   'printsiteals 1',
                                   ]))
    else:
        paramFile.write('\n'.join(['Mode WMREF',
                                   'refspecies hg19',
                                   'TREE (hg19:1);',
                                   'wmdiff 0.05',
                                   'EMprior 1',
                                   'priordiff 0.01',
                                   'markovorderBG 0',
                                   'bgprior 0.999',
                                   'bg A %0.5f' % args.background['A'],
                                   'bg C %0.5f' % args.background['C'],
                                   'bg G %0.5f' % args.background['G'],
                                   'bg T %0.5f' % args.background['T'],
                                   'restrictparses 0',
                                   'sitefile %s' % os.path.join(args.output_dir, 'sites'),
                                   'priorfile %s' % os.path.join(args.output_dir, 'priors'),
                                   'minposterior %f' % args.min_post,
                                   'printsiteals 1',
                                   ]))
    paramFile.flush()
    return 0


def copyWM(args, TF, program_dir=''):
    motifLength = len([line for line in open(args.WM)
                       if re.search(r'^\d+', line)])
    motevoParamFile(os.path.join(args.output_dir, 'motevo_params'), motifLength, TF, args)
    runMotevo(args.fasta_file,
              'motevo_params',
              os.path.basename(args.WM), args.output_dir, program_dir)
    return 0


def process_motevo_output(args):
    site_file = os.path.join(args.output_dir, 'sites')
    alg = []
    with open(os.path.join(args.output_dir, 'alg_1'), 'w') as outf:
        with open(site_file) as inf:
            while True:
                line1 = inf.readline()
                if not line1:
                    break
                line2 = inf.readline()
                outf.write('%s\t%s\n' % (line2.split()[0], line1.split()[2]))
                alg.append( (line2.split()[0], float(line1.split()[2])) )
    return os.path.join(args.output_dir, 'alg_1'), alg


def pairwise_frequencies(sequences, length):
    convert = {'A':0, 'C':1, 'G':2, 'T':3}
    dinucleotide_matrix = np.zeros(length*length*16).reshape(length, length, 16)
    for s in sequences:
        for p in combinations(range(length), 2):
            dinucleotide_matrix[p[0]][p[1]][ convert[s[0][p[0]].upper()]*4 + convert[s[0][p[1]].upper()] ] += s[1]
    return dinucleotide_matrix


def generate_DWT(alg_file, TF):
    args = arguments()
    # test if the input file contains two columns or only one
    sequences = []
    f = open(alg_file)
    if len(f.readline().split('\t')) == 2:
        f.seek(0)
        sequences = [(s.split()[0], float(s.split()[1])) for s in f]
    else:
        f.seek(0)
        sequences = [(s.split()[0], 1.0) for s in f]

    length = len(sequences[0][0])
    dinucleotide_matrix = pairwise_frequencies(sequences, length)
    header = [''.join(i) for i in product('ACGT',repeat=2)]

    dwt_file = os.path.join(args.output_dir, os.path.basename(alg_file).replace('alg', 'dwt'))
    with open(dwt_file, 'w') as outf:
        if TF is None:
            outf.write('NA\tTF\n')
        else:
            outf.write('NA\t%s\n' % (TF))

        outf.write('PO1\tPO2\t%s\n' % ('\t'.join(header)))
        for p in combinations(range(length), 2):
            outf.write('%d\t%d\t%s\n' % (p[0]+1, p[1]+1, '\t'.join(map(str,
                                                                 np.round(dinucleotide_matrix[p[0]][p[1]][:],3) ) ) ))
    return dwt_file


def background_model(infile):
    """
    As input, this method takes a filename which contains the sequences in FASTA format
    and returns the sequences' base frequencies in a form of a dictionary.
    """
    test = re.compile ("^>>", re.IGNORECASE)
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
    for b in freq.keys(): freq[b] = freq[b]/normalization
    return freq


def DWT_param_file(args, res_filename):
    if args.with_background:
        str = '\n'.join(['#background',
                         'A\t%0.5f' % args.background['A'],
                         'C\t%0.5f' % args.background['C'],
                         'G\t%0.5f' % args.background['G'],
                         'T\t%0.5f' % args.background['T'],
                         '#minimum_score',
                         '-80.0000',
                         '#output_file',
                         '%s' % res_filename
                         ])
    else:
        str = '\n'.join(['#background',
                         'A\t0.25',
                         'C\t0.25',
                         'G\t0.25',
                         'T\t0.25',
                         '#minimum_score',
                         '-80.0000',
                         '#output_file',
                         '%s' % res_filename
                         ])

    fname = os.path.join(args.output_dir, 'DWT_param_file')
    try:
        f = open(fname, 'w')
        f.write(str)
        f.close()
    except IOError:
        print 'Error in making %s' % f_name
        return False
    return fname


def run_DWT_model(param_file, dwt_file, args, program_dir):
    prog = os.path.join(program_dir, '../DWT_model/Source/DWT_TFBS_prediction')
    cmd = '%s %s %s %s' % (prog, dwt_file, args.fasta_file, param_file)
    if run(cmd):
         sys.stderr.write ( "\nError in running the DWT model\n" )
         print "command: %s" % cmd
         print 'Program halts!\n'
         exit()
    return 0

def pr(s_i, s_0):
    return np.exp(s_i - s_0) / (1 + np.exp(s_i - s_0))

def find_cutoff(S):
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

def freq_matrix_4D(alg):
    """returns a matrix of the letter frequences"""
    matrix = {}
    for seq in alg:
        for pos in combinations(xrange(len(alg[0][0])), 2):
            matrix.setdefault(pos, {})
            for pairs in product('ACGT', repeat=2):  matrix[pos].setdefault(pairs, 0.)
            try:
                matrix[pos][(upper(seq[0][pos[0]]), upper(seq[0][pos[1]]))] += seq[1]
            except KeyError:
                continue
    return matrix


def difference(alg1, alg2=None):
    if alg2 is None:
        return 10000000.   # just a big number to stay above the epsilon in iterate() method

    if len(alg1[0][0]) != len(alg2[0][0]):
        raise Exception('UnequalLengths')
    (matrix1, matrix2) = map(freq_matrix_4D, [alg1, alg2])
    sum = 0.
    for pos in combinations(xrange(max([len(alg1[0][0]), len(alg2[0][0])])), 2):
        for pairs in product('ACGT', repeat=2):
            sum += 2*math.fabs(matrix1[pos][pairs] - matrix2[pos][pairs]) # 2 * |N_ij^AB - M_ij^AB|
            if (matrix1[pos][pairs] + matrix2[pos][pairs]):
                sum /= (matrix1[pos][pairs] + matrix2[pos][pairs])
    return sum


def make_alignment(dwt_res_file, args, s_0=None, iteration=0):
    if s_0 is None:
        scores = [float(a_line.split()[-1]) for a_line in open(dwt_res_file)]
        s_0 = find_cutoff(scores)
        del(scores)
    # print "s0 is: ", s_0
    alg = []
    for line in open(dwt_res_file):
        if pr(float(line.split()[-1]), s_0) > args.min_post:
            alg.append((line.split()[-2], pr(float(line.split()[-1]), s_0)))
    alg_file = os.path.join(args.output_dir, 'alg_%d' % iteration)
    with open(alg_file, 'w') as outf:
        for rec in alg:
            outf.write('%s\t%f\n' % (rec[0], rec[1]))
    # print "Number of sequences in the alignment file is %d" % len(alg)
    return alg, s_0, alg_file


def clean_up_directory(args, last_iteration, TF):
    curr_dir = os.getcwd()
    os.chdir(args.output_dir)
    run('mkdir intermediate_results')
    run('mv * intermediate_results/.')
    run('cp intermediate_results/dwt_%d %s.dwt' % (last_iteration, TF))
    run('cp intermediate_results/alg_%d %s.alg' % (last_iteration, TF))
    os.chdir(curr_dir)
    return 0


def calculate_PS_posterior(args, TF, program_dir=''):
    prog = os.path.join(program_dir, '../Positional_Dependency_Posterior/Source/positional_dependency_posterior')
    cmd = ' '.join([
        prog,
        os.path.join(args.output_dir,'%s.dwt' % TF),
        '>',
        os.path.join(args.output_dir,'%s.post' % TF),
    ])
    if run(cmd):
         sys.stderr.write ( "\nError in calculating positional dependencies\n" )
         print "command: %s" % cmd
         print 'Program halts!\n'
         exit()
    return 0


def main():
    args =arguments()
    mkdir(args.output_dir)
    program_dir = os.path.dirname(sys.argv[0])
    if not args.TF:
        TF = os.path.basename(args.WM).split('.')[0]
    else:
        TF = args.TF

    if args.with_background:
        args.background = background_model(args.fasta_file)
    else:
        args.background = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    iteration = 1
    copyWM(args, TF, program_dir)
    alg_file, alg_1 = process_motevo_output(args)
    print "number of identified TFBS in iteration 1 is %d" % len(alg_1)
    print "*****************************"
    dwt_file = generate_DWT(alg_file, TF)
    dwt_results = os.path.join(args.output_dir, 'DWT_results')
    dwt_param_file = DWT_param_file(args, dwt_results)
    diff = 100
    while diff > 0.01:
        iteration += 1
        run_DWT_model(dwt_param_file, dwt_file, args, program_dir)
        alg_2, s0, alg_file2 = make_alignment(dwt_results, args, iteration=iteration)
        dwt_file = generate_DWT(alg_file2, TF)
        diff = difference(alg_1, alg_2)
        print "difference at round %i is %f" % (iteration, diff)
        print "number of identified TFBS in this iteration is %d" % len(alg_2)
        print "*****************************"
        alg_1 = alg_2
        alg_2 = None
    clean_up_directory(args, 12, TF)
    calculate_PS_posterior(args, TF, program_dir)

if __name__ == '__main__':
        main()