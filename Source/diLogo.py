
## Constants 
left_offset = 2.
low_offset = 7.
width = 1.5
distance = 2.
height = 5.
radius = 0.75
pseudo_count = .5
bp_index = {0:'T', 1:'G', 2:'C', 3:'A'}


def arguments():
    import argparse
    parser = argparse.ArgumentParser(description="""
    For generating graphical representation of the binding specificity under
    the DWT model, this code takes as input a DWT flat file and a file containing
    positional dependency posteriors between all pairs of positions. In return,
    it generates a diLogo in PDF format, located either at the current directory
    or a directory that is given.
    """)

    parser.add_argument('-i', dest='dwt_file', action='store', required=True, \
                        help="""A DWT flat file that encodes the binding specificity of a
                        TF under the DWT model. """)
    parser.add_argument('-p', dest='posterior_file', action='store', required=True, \
                        help="""The posteriors between all pair of positions.
                        A C++ program, located in 'Positional_Dependency_Posterior' can be
                        used to calculate the posteriors of dependencies.""")
    parser.add_argument('-c', dest='cutoff', action='store', required=False, default=0., type=float,  \
                        help="""Optional cutoff over the posteriors of dependency. This
                        argument is used to adjust the number of dependencies in the final
                        diLogo. By default cutoff is set to 0.""")
    parser.add_argument('-o',  dest='out_dir', action='store',
                        type=str, required=False,
                        help="""Optional argument that gives the output directory in which
                        the resulting PDF file will be copied in. The directory must exist.
                        If not given, the current working directory is used.""")
    args = parser.parse_args()
    return args


def make_spanning_tree(sorted_posteriors):
    spanning_tree = [sorted_posteriors[0][0]]
    n = int((sqrt(1+8*len(sorted_posteriors)) + 1)/2)
    for i in xrange(n-2):
        for e in sorted_posteriors:
            if no_loop(spanning_tree,e[0]):
                spanning_tree.append(e[0])
                break
    return spanning_tree


def insert(bp, x, y, height):    
    c.insert(bitmap.bitmap(x, y, nucleotide[bp], height=height, width=width, compressmode=None))        


def no_loop(tree, e):
    tree_nodes = set(reduce(add, tree))
    if len(tree_nodes.intersection(set(e))) == 1:
        return True
    return False

    
def prune_tree(tree, cutoff, posteriors):
    """
    removes edges with the posterior lower than the cutoff
    """
    new_tree = []
    for e in tree:
        try:
            if posteriors[e] > cutoff:
                new_tree.append(e)
        except KeyError:
            if posteriors[e[::-1]] > cutoff:
                new_tree.append(e)
    return new_tree


def change_directions(tree):
    """
    The function arrange the directions of the edges in a way that every node
    gain in-degree of one    
    """
    tmp = [] # holds the nodes that have edges pointing to
    new_tree = []
    for e in tree:
        try:
            if tmp.index(e[1])>=0:
                new_tree.append(e[::-1])
                tmp.append(e[0])
        except ValueError:
            new_tree.append(e)
            tmp.append(e[1])
    return new_tree


def directions(tree, root):
    """
    The function arrange the directions of the edges in a way that every node
    gain in-degree of one    
    """
    tmp = [] # holds the nodes that have edges pointing to    
    new_tree = []
    for e in tree:
        if e[0] == root:
            new_tree.append((root,e[1]))            
            tmp.append(e[1])
            tree.remove(e)
        elif e[1] == root:
            new_tree.append((root, e[0]))
            tmp.append(e[0])
            tree.remove(e)
            
    for e in tree:
        try:            
            if tmp.index(e[1])>=0:
                new_tree.append(e[::-1])
                tmp.append(e[0])
        except ValueError:
            new_tree.append(e)
            tmp.append(e[1])
    return new_tree


def select_root(posteriors):
    v = []
    for i in xrange(n):
        v.append(0.)
        for j in xrange(n):
            if i < j:
                v[i] += posteriors[(i,j)]
            elif i > j:
                v[i] += posteriors[(j,i)]
    return v.index(max(v))


def freq_matrix_4D(dwt_file):
    """returns a matrix of the letter frequences"""
    matrix = {}       
    for line in open(dwt_file):
        if re.search('PO', line):
            letters = [(i[0], i[1]) for i in line.split()[2:]]
        if not (re.search('PO', line) or re.search('NA', line) or line.startswith('\\')):
            tmp = line.split('\t')
            pos = (int(tmp[0])-1, int(tmp[1])-1)
            matrix.setdefault(pos, {})
            for i in xrange(16):
                matrix[pos][letters[i]] = float(tmp[i+2])
    return matrix


def probabilities_matrix_4D(count_matrix, n):
    matrix = {}
    for pos in count_matrix.keys():
        for bp in product('ACGT', repeat=2):
            matrix.setdefault(pos, {})
            matrix[pos][bp] = count_matrix[pos][bp] / float(n)
    return matrix


def probabilities_matrix_2D(count_matrix, n):
    matrix = {}
    for pos in count_matrix.keys():
        for bp in 'ACGT':
            matrix.setdefault(pos, {})
            matrix[pos][bp] = count_matrix[pos][bp] / float(n)
    return matrix
    
    
def freq_matrix_2D(matrix_4D, n):
    """returns a matrix of the letter frequences"""
    matrix = {}
    for i in xrange(n):
        matrix.setdefault(i, {})
        if i < n-1:
            for bp in 'ACGT':
                matrix[i][bp] = matrix_4D[(i, i+1)][(bp, 'A')] + matrix_4D[(i, i+1)][(bp, 'C')] +  matrix_4D[(i, i+1)][(bp, 'G')] + matrix_4D[(i, i+1)][(bp, 'T')]
        else:
            for bp in 'ACGT':
                matrix[i][bp] = matrix_4D[(i-1, i)][('A', bp)] + matrix_4D[(i-1, i)][('C', bp)] + matrix_4D[(i-1, i)][('G', bp)] + matrix_4D[(i-1, i)][('T', bp)]

    return matrix


def draw(pos, bp, freq, N, sum):
    entropy = -1*reduce(lambda x,y: x+y,
                        map(lambda f: f*log(f, 2),freq.values()))
    heights = {}
    # correction = 3./(2*log(2)*N)
    correction = 0.    
    for n in 'ACGT':
        heights[n] = freq[n]*(2-(entropy+correction))
    sorted_heights = sorted(heights.items(), key=itemgetter(1))
    y_shift = height/20
    for i in sorted_heights:
        if i[1] <=  (pseudo_count/sum)*(2-(entropy+correction)):        
        # if i[1] < 0.05:
            continue
        insert(i[0], pos_x[pos[1]]-(width/2), pos_y[bp]+y_shift, (height/2)*i[1])
        y_shift += (height/2)*i[1] + (height/50) 

def TF_name(dwt_file):
    TF = 'TF'
    for line in open(dwt_file):
        if re.search('NA', line):
            TF = line.split('\t')[-1].rstrip()
    return TF


def draw_WM(x, freq, N):
    # entropy = -1*reduce(lambda x,y: x+y,
    #                     map(lambda f: f*log(f, 2),freq.values()))

    entropy = 0
    for i in freq.values():
        try:
            entropy -= i*log(i, 2)
        except ValueError:
            continue
    heights = {}
    # correction = 3./(2*log(2)*N)
    correction = 0.
    for n in 'ACGT':
        heights[n] = freq[n]*(2 - (entropy+correction))
    sorted_heights = sorted(heights.items(), key=itemgetter(1))
    # y_shift = height/25
    y_shift = 0
    for i in sorted_heights:
        if i[1] <=  (pseudo_count/sum)*(2-(entropy+correction)):
            continue
        insert(i[0], pos_x[x]-(width/2), pos_y['A']+ height+y_shift+height/10, (height/2)*i[1])
        y_shift += ((height/2)*i[1] + height/50)
    return True

if __name__ == '__main__':
    import sys, re, os
    from itertools import *
    from math import *
    from pyx import *
    from operator import itemgetter, add
    from string import upper

    code_dir = os.path.dirname(sys.argv[0])

    args = arguments()
    posteriors = dict([((int(i.split()[0]), int(i.split()[1])), float(i.split()[2]))
                       for i in open(args.posterior_file ) 
                       if not re.search(r'alpha_exponent', i)])
    sorted_posteriors = sorted(posteriors.items(), key=itemgetter(1), reverse=True)
    n = int((sqrt(1+8*len(sorted_posteriors)) + 1)/2)
    spanning_tree = make_spanning_tree(sorted_posteriors)
    TF = TF_name(args.dwt_file)
    text.set(mode="latex")
    text.preamble(r"\renewcommand{\familydefault}{\sfdefault}")
    c = canvas.canvas()
    pos_y = {}
    for i in xrange(4):
        c.stroke(path.rect(left_offset, low_offset+i*(height+(height/25)), n*distance + (distance/5), height+(height/25)),
                 [style.linewidth.Thick])
        pos_y[bp_index[i]] = low_offset+i*(height+(height/25))

    c.stroke(path.line(left_offset-0.3, pos_y['A']+height+(height/25)-0.04, left_offset-0.3, pos_y['A']+2*(height+(height/25))+0.04), \
             [style.linewidth(0.1)])
    c.stroke(path.line(left_offset-0.3, pos_y['A']+2*(height+(height/25)),
                       left_offset-0.8, pos_y['A']+2*(height+(height/25))), [style.linewidth(0.1)])
    c.stroke(path.line(left_offset-0.3, pos_y['A']+(height+(height/25)),
                       left_offset-0.8, pos_y['A']+(height+(height/25))), [style.linewidth(0.1)])
    c.stroke(path.line(left_offset-0.3, pos_y['A']+1.5*(height+(height/25)),
                       left_offset-0.8, pos_y['A']+1.5*(height+(height/25))), [style.linewidth(0.1)])
    c.text(left_offset-width-(width/3), pos_y['A']+1.5*(height+(height/25)), r"bits", [text.halign.boxcenter, text.size(5), trafo.rotate(90)])
    c.text(left_offset-1.7, pos_y['A']+2*(height+(height/25))-.1, r"2.0", [text.halign.boxcenter, text.size(3)])
    c.text(left_offset-1.7, pos_y['A']+(height+(height/25))-.1, r"0.0", [text.halign.boxcenter, text.size(3)])

    c.text((left_offset + n*distance + (distance/5))/2, pos_y['A']+2.2*(height+(height/25)), r"\textbf{%s}" % TF,
       [text.size(5), text.halign.boxcenter])
         
    pos_x = {}

    for i in xrange(n):
        c.stroke(path.circle(left_offset+i*distance+radius*2, low_offset-radius-(radius/2), radius),
                 [color.cmyk.RoyalBlue,
                  deco.filled([color.cmyk.CornflowerBlue])])
        c.text(left_offset+i*distance+radius*2, low_offset-2*radius + radius/5, r"%d" % (i+1), [text.halign.boxcenter,
                                                                    text.size(3)])                                                           
        pos_x[i] = left_offset+i*distance + radius*2

        nucleotide = dict([(bp, bitmap.jpegimage(os.path.join(code_dir, "nucleotides_images/%s.jpg" % bp)))
                           for bp in 'ACGT'])

        for i,bp in zip(range(4), 'TGCA'):    
            insert(bp, left_offset-width-(width/5), low_offset+i*(height+(height/25))+(height/20), height-(height/4))
            
    root = select_root(posteriors)
    print "root is: %d" % (root+1)
    tree = directions(spanning_tree, root)
    tree = prune_tree(tree, args.cutoff, posteriors)
    print tree

    matrix_4D = freq_matrix_4D(args.dwt_file)
    matrix_2D = freq_matrix_2D(matrix_4D, n)
    total_number_of_sequences = matrix_2D[0]['A'] + matrix_2D[0]['C'] +matrix_2D[0]['G'] + matrix_2D[0]['T']
    
    freq_matrix = probabilities_matrix_4D(matrix_4D, total_number_of_sequences)
    for x in tree:
        if x[0]>x[1]:
            pos = x[::-1]
            for n in 'ACGT':
                sum = 4*pseudo_count
                for b in 'ACGT':
                    sum += matrix_4D[pos][(b,n)]
                freq = {}
                for b in 'ACGT':
                    freq[b] = (matrix_4D[pos][(b,n)] + pseudo_count)/sum
                    # freq[b] = (freq_matrix[pos][(b,n)])/sum                
                for i in freq.keys():
                    draw(pos[::-1], n, freq, total_number_of_sequences, sum)
        else:
            pos = x
            for n in 'ACGT':
                sum = 4*pseudo_count
                for b in 'ACGT':
                    sum += matrix_4D[pos][(n,b)]
                freq = {}
                for b in 'ACGT':
                    freq[b] = (matrix_4D[pos][(n,b)] + pseudo_count)/sum
                    # freq[b] = (freq_matrix[pos][(n,b)])/sum                               
                for i in freq.keys():                
                    draw(pos, n, freq, total_number_of_sequences, sum)

    freq_matrix = probabilities_matrix_2D(matrix_2D, total_number_of_sequences)
    n = int((sqrt(1+8*len(sorted_posteriors)) + 1)/2)
    for pos in xrange(n):
        sum = 4*pseudo_count
        for b in 'ACGT':
            sum += matrix_2D[pos][b]
        freq = {}    
        for b in 'ACGT':
            freq[b] = (matrix_2D[pos][b] + pseudo_count)/sum
        draw_WM(pos, freq, total_number_of_sequences)

    # put the edges
    curve_smooth = radius/2
    from random import random
    for e in sorted(tree):
        c.stroke(path.curve(pos_x[e[0]]+random()/8, low_offset-2*radius-(radius/2),
                            pos_x[e[0]], low_offset-2-curve_smooth,
                            pos_x[e[1]], low_offset-2-curve_smooth,
                            pos_x[e[1]]-random()/6, low_offset-2*radius-(radius/2)),                        
                 [color.cmyk.RoyalBlue,              
                  style.linewidth(0.07),
                  deco.earrow.Large()])
        curve_smooth += radius/3

    # heat map

    for i in xrange(0,n):
        for j in xrange(i+1,n):
            c.stroke(path.rect(pos_x[j]-(width/2)-(width/6), distance*(-1*(i+2)) +  low_offset/2, distance, distance),
                     [deco.filled([color.rgb.red, color.transparency(1.0 - posteriors[(i,j)])]), color.cmyk.RedOrange])
    for i in xrange(0,n-1):    
        c.text(pos_x[0] + i*distance, distance*(-1*(i+2))+low_offset/2 + distance/2.3, r"%d" % (i+1),
               [text.size(4), text.halign.boxcenter])  # left side
    bitmap.jpegimage(os.path.join(code_dir, "nucleotides_images/%s.jpg" % bp))
    c.insert(bitmap.bitmap(pos_x[0] - 2.2*distance, distance*(-1*n) + low_offset/2 - distance/5,
                           bitmap.jpegimage(os.path.join(code_dir, "nucleotides_images/key.jpg")),
                            height=distance*(n-1) + distance/2, width=distance*1.1, compressmode=None, dctquality=100))
    if not args.out_dir:
        c.writePDFfile(TF)
    else:
        c.writePDFfile(os.path.join(args.out_dir, TF))

