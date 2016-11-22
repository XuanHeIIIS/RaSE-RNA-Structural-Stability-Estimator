#!/usr/bin/env python


"""RaSE - RNA structurAl Stability Estimate.

Compute stability.

Version: 1.0
Author: Fabrizio Costa [costa@informatik.uni-freiburg.de]

Usage:
  RaSE -i <sequence>
              [-k N] [-c N] [-n N] [-w N] [-b N] [-p N] [-r N] [-e N]
              [-l] [-t] [--draw]
              [--verbose]
  RaSE (-h | --help)
  RaSE --version

Options:
  -i <sequence>                     Specify input sequence.
  -k N                              Specify number of maximally unstable
                                    nucleotides to mark [default: 5].
  -c N                              Complexity of features [default: 3].
  -n N                              Size feature pseudo identifiers
                                    [default: 15].
  -w N                              Window size [default: 150]
  -b N                              Max number of spanning bases [default: 130]
  -p N                              Average probability cutoff [default: 0.1]
  -r N                              Hard threshold [default: 0.5]
  -e N                              Max num edges [default: 2]
  -l                                Flag to activate no lonely base pairs mode.
  -t                                Flag to activate no nesting mode.
  --draw                            Output drawing with standard name out.pdf.
  -h --help                         Show this screen.
  --version                         Show version.
  --verbose                         Print more text.


"""
from docopt import docopt
from toolz import curry, compose
from toolz.sandbox.core import unzip
from sklearn.metrics.pairwise import pairwise_kernels
from eden.converter.rna.rnaplfold import rnaplfold_to_eden
from eden.graph import Vectorizer
from eden.util.display import draw_graph
from eden.util import configure_logging, serialize_dict

import logging
logger = logging.getLogger(__name__)


def make_fold(window_size=150,
              max_bp_span=130,
              avg_bp_prob_cutoff=0.1,
              hard_threshold=0.5,
              max_num_edges=2,
              no_lonely_bps=True,
              nesting=True):
    fold = curry(rnaplfold_to_eden)(window_size=window_size,
                                    max_bp_span=max_bp_span,
                                    avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                                    hard_threshold=hard_threshold,
                                    max_num_edges=max_num_edges,
                                    no_lonely_bps=no_lonely_bps,
                                    nesting=nesting)
    return fold


def make_fold_vectorize(complexity=3,
                        nbits=15,
                        fold=None):
    vec = Vectorizer(complexity=complexity, nbits=nbits)
    vectorize = curry(lambda vec, graphs: vec.transform(graphs))(vec)

    fold_vectorize = compose(vectorize, fold)
    return fold_vectorize


def make_string(gen):
    return ''.join(list(gen))


def make_variations(seq, index=0, alphabet='ACGU'):
    for nt in alphabet:
        snip = list(seq)
        if snip[index] != nt:
            snip[index] = nt
            snip = make_string(snip)
            yield nt, snip


def compute_stability(seq, alphabet='ACGU', fold_vectorize=None):
    for index in range(len(seq)):
        cmake_variations = curry(make_variations)(index=index,
                                                  alphabet=alphabet)
        make_vecs = compose(fold_vectorize, cmake_variations)
        variation_vecs = make_vecs(seq)

        vec = fold_vectorize([('', seq)])

        sims = pairwise_kernels(vec, variation_vecs, metric='cosine')
        sims = sims.reshape((len(alphabet) - 1, 1))
        sims_variations = zip(sims, cmake_variations(seq))

        score, snip = min([(sim[0], variation[0])
                           for sim, variation in sims_variations])

        yield score, snip


def stability(seq, alphabet='ACGU', fold_vectorize=None):
    scores, snips = unzip(
        compute_stability(seq,
                          alphabet=alphabet,
                          fold_vectorize=fold_vectorize))
    scores = list(scores)
    snips = make_string(snips)
    return snips, scores


def serialize(seq, snips, scores, k=5):
    tuples = sorted(
        [(score, i, nt, snip)
         for i, (snip, score, nt) in enumerate(zip(snips, scores, seq))])
    score_th = tuples[k][0]
    for i, (snip, score, nt) in enumerate(zip(snips, scores, seq)):
        if score <= score_th:
            mark = '*'
        else:
            mark = ''
        yield '%3d %s %s %.2f %s' % (i, nt, snip, score, mark)


def rna_structural_stability_estimate(seq,
                                      alphabet='ACGU',
                                      k=5,
                                      complexity=3,
                                      nbits=15,
                                      fold=None):
    fold_vectorize = make_fold_vectorize(complexity=complexity,
                                         nbits=nbits,
                                         fold=fold)
    snips, scores = stability(seq, alphabet='ACGU',
                              fold_vectorize=fold_vectorize)
    return serialize(seq, snips, scores, k=k)


def update_graph(graph, scores, snips):
    for u in graph.nodes():
        pos = graph.node[u]['position']
        graph.node[u]['stability'] = scores[pos]
        graph.node[u]['snip'] = snips[pos]
    return graph


def draw(seq,
         file_name=None,
         fold=None,
         complexity=3,
         nbits=15):
    fold_vectorize = make_fold_vectorize(complexity=complexity,
                                         nbits=nbits,
                                         fold=fold)
    snips, scores = stability(seq, fold_vectorize=fold_vectorize)

    seqs = [('', seq)]
    graphs = fold(seqs)
    graph = graphs.next()
    graph = update_graph(graph, scores, snips)

    opts = {'size': 14, 'vertex_border': False, 'vertex_size': 400,
            'edge_label': None, 'font_size': 9, 'vertex_alpha': 0.6,
            'invert_colormap': True, 'secondary_vertex_label': 'snip',
            'vertex_color': 'stability', 'colormap': 'Blues',
            'ignore_for_layout': 'nesting'}
    draw_graph(graph, file_name=file_name, **opts)


def main(args):
    """Main."""
    # setup variables
    seq = args['-i']
    k = int(args['-k'])
    complexity = int(args['-c'])
    nbits = int(args['-n'])
    window_size = int(args['-w'])
    max_bp_span = int(args['-b'])
    avg_bp_prob_cutoff = float(args['-p'])
    hard_threshold = float(args['-r'])
    max_num_edges = int(args['-e'])
    no_lonely_bps = args['-l']
    no_nesting = args['-t']
    if no_nesting is True:
        nesting = False
    else:
        nesting = True
    # logger
    if args['--verbose']:
        verbosity = 2
    else:
        verbosity = 1
    configure_logging(logger,
                      verbosity=verbosity,
                      filename='log')
    logger.debug(serialize_dict(args))

    fold = make_fold(window_size=window_size,
                     max_bp_span=max_bp_span,
                     avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                     hard_threshold=hard_threshold,
                     max_num_edges=max_num_edges,
                     no_lonely_bps=no_lonely_bps,
                     nesting=nesting)
    rase = curry(rna_structural_stability_estimate)(k=k,
                                                    complexity=complexity,
                                                    nbits=nbits,
                                                    fold=fold)
    for line in rase(seq):
        print(line)

    if args['--draw']:
        draw(seq,
             fold=fold,
             complexity=complexity,
             nbits=nbits,
             file_name='out.pdf')

if __name__ == '__main__':
    args = docopt(__doc__, version='RaSE v1.0')
    main(args)
