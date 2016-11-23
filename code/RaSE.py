#!/usr/bin/env python


"""RaSE - RNA structurAl Stability Estimate.

Compute stability.

Version: 1.0
Author: Fabrizio Costa [costa@informatik.uni-freiburg.de]

Usage:
  RaSE [-i <sequence>]
       [-k N] [-c N, --complexity=N] [-n N, --nbits=N] [-w N, --window_size=N]
       [-b N, --max_bp_span=N] [-p N, --avg_bp_prob_cutoff=N]
       [-r N, --hard_threshold=N] [-e N, --max_num_edges=N]
       [-l, --no_lonely_bps] [-t, --no_nesting] [--draw]
       [--verbose]
  RaSE (-h | --help)
  RaSE --version

Options:
  -i <sequence>                     Specify input sequence [default: stdin].
  -k N                              Specify number of maximally unstable
                                    nucleotides to mark [default: 5].
  -c N, --complexity=N              Complexity of features [default: 3].
  -n N, --nbits=N                   Num bits to represent all possible feature
                                    pseudo identifiers [default: 15].
  -w N, --window_size=N             Window size [default: 150]
  -b N, --max_bp_span=N             Max number of spanning bases [default: 130]
  -p N, --avg_bp_prob_cutoff=N      Average probability cutoff [default: 0.1]
  -r N, --hard_threshold=N          Hard threshold [default: 0.5]
  -e N, --max_num_edges=N           Max num edges [default: 2]
  -l, --no_lonely_bps               Flag to activate no lonely base pairs mode.
  -t, --no_nesting                  Flag to activate no nesting mode.
  --draw                            Output drawing with standard name out.pdf.
  -h --help                         Show this screen.
  --version                         Show version.
  --verbose                         Print more text.


"""
import sys
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
    """Curry parameters in the folding algorithm.

     Parameters
    ----------
    max_num_edges : int (default 2)
        Is the maximal number of base pair bonds allowed per nucleotide.

    window_size : int (default 150)
        Is the size of the window.

    max_bp_span : int (default 130)
        Is the maximum number of bases between to bases that pair.

    avg_bp_prob_cutoff : float (default 0.1)
        Is the threshold value under which the edge is not materialized.

    no_lonely_bps : bool (default True)
        If True no lonely base pairs are allowed.

    nesting : bool (default True)
        If True the edge type is 'nesting'.

    hard_threshold : float (default 0.5)
        If the edges with avg probability is greater than hard_threshold
        then they are not of nesting type.
    """
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
    """Curry parameters in vectorizer."""
    vec = Vectorizer(complexity=complexity, nbits=nbits)
    vectorize = curry(lambda vec, graphs: vec.transform(graphs))(vec)

    fold_vectorize = compose(vectorize, fold)
    return fold_vectorize


def _make_string(gen):
    return ''.join(list(gen))


def _make_variations(seq, index=0, alphabet='ACGU'):
    for nt in alphabet:
        snip = list(seq)
        if snip[index] != nt:
            snip[index] = nt
            snip = _make_string(snip)
            yield nt, snip


def compute_stability(seq, alphabet='ACGU', fold_vectorize=None):
    """Compute the structural effects of single nt change.

    Specifically, compute the least similarity of the structure obtained
    by replacing each nucleotide with all possible alternatives.
    """
    # TODO: parallelize indices
    for index in range(len(seq)):
        _cmake_variations = curry(_make_variations)(index=index,
                                                    alphabet=alphabet)
        make_vecs = compose(fold_vectorize, _cmake_variations)
        variation_vecs = make_vecs(seq)

        vec = fold_vectorize([('', seq)])

        sims = pairwise_kernels(vec, variation_vecs, metric='cosine')
        sims = sims.reshape((len(alphabet) - 1, 1))
        sims_variations = zip(sims, _cmake_variations(seq))

        score, snip = min([(sim[0], variation[0])
                           for sim, variation in sims_variations])

        yield score, snip


def stability(seq, alphabet='ACGU', fold_vectorize=None):
    """Compute the structural effects of single nt change.

    Specifically, compute the least similarity of the structure obtained
    by replacing each nucleotide with all possible alternatives.
    This function wraps compute_stability and post processes its output.
    """
    scores, snips = unzip(
        compute_stability(seq,
                          alphabet=alphabet,
                          fold_vectorize=fold_vectorize))
    scores = list(scores)
    snips = _make_string(snips)
    return snips, scores


def serialize(seq, snips, scores, k=5):
    """Pretty print of the snip information and the realtive scores."""
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
    """Wrapper for the computation of the structural effects of nt change.

    This function pretty prints the computed results.
    """
    fold_vectorize = make_fold_vectorize(complexity=complexity,
                                         nbits=nbits,
                                         fold=fold)
    snips, scores = stability(seq, alphabet='ACGU',
                              fold_vectorize=fold_vectorize)
    return serialize(seq, snips, scores, k=k)


def _update_graph(graph, scores, snips):
    for u in graph.nodes():
        pos = graph.node[u]['position']
        graph.node[u]['stability'] = scores[pos]
        graph.node[u]['snip'] = snips[pos]
    return graph


def _draw(seq, file_name=None, snips=None, scores=None, fold=None):
    seqs = [('', seq)]
    graphs = fold(seqs)
    graph = graphs.next()
    graph = _update_graph(graph, scores, snips)

    opts = {'size': 14, 'vertex_border': False, 'vertex_size': 400,
            'edge_label': None, 'font_size': 9, 'vertex_alpha': 0.6,
            'invert_colormap': True, 'secondary_vertex_label': 'snip',
            'vertex_color': 'stability', 'colormap': 'Blues',
            'ignore_for_layout': 'nesting', 'layout': 'KK'}
    draw_graph(graph, file_name=file_name, **opts)


def draw(seq,
         file_name=None,
         fold=None,
         complexity=3,
         nbits=15):
    """Graphical representation of the RNA folded structure.

    Nodes encode the original nt and the most de-stabilizing alternative.
    """
    fold_vectorize = make_fold_vectorize(complexity=complexity,
                                         nbits=nbits,
                                         fold=fold)
    snips, scores = stability(seq, fold_vectorize=fold_vectorize)
    _draw(seq, file_name=file_name, snips=snips, scores=scores, fold=fold)


def main(args):
    """Main."""
    # read variables
    # if no -i is given then read from stdin
    seq = args['-i']
    seq = (sys.stdin.readline().strip() if args['-i'] == 'stdin' else seq)
    k = int(args['-k'])
    complexity = int(args['--complexity'][0])
    nbits = int(args['--nbits'][0])
    window_size = int(args['--window_size'][0])
    window_size = min(len(seq), window_size)
    max_bp_span = int(args['--max_bp_span'][0])
    max_bp_span = min(len(seq), max_bp_span)
    avg_bp_prob_cutoff = float(args['--avg_bp_prob_cutoff'][0])
    hard_threshold = float(args['--hard_threshold'][0])
    max_num_edges = int(args['--max_num_edges'][0])
    no_lonely_bps = args['--no_lonely_bps']
    no_nesting = args['--no_nesting']
    if no_nesting is True:
        nesting = False
    else:
        nesting = True
    # setup logger
    if args['--verbose']:
        verbosity = 2
    else:
        verbosity = 1
    configure_logging(logger,
                      verbosity=verbosity,
                      filename='log')
    logger.debug(serialize_dict(args))

    # setup folding algorithm
    fold = make_fold(window_size=window_size,
                     max_bp_span=max_bp_span,
                     avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                     hard_threshold=hard_threshold,
                     max_num_edges=max_num_edges,
                     no_lonely_bps=no_lonely_bps,
                     nesting=nesting)
    # setup vectorizer
    fold_vectorize = make_fold_vectorize(complexity=complexity,
                                         nbits=nbits,
                                         fold=fold)
    # compute the most de-stabilizing single nt change per position
    # and the corresponding similarity of the resulting structure
    # w.r.t. the original
    snips, scores = stability(seq,
                              alphabet='ACGU',
                              fold_vectorize=fold_vectorize)
    # print: nt pos, original nt, most de-stabilizing nt, similarity score
    for line in serialize(seq, snips, scores, k=k):
        print(line)

    # if drawing is required use the folding algorithm to compute the graph
    if args['--draw']:
        _draw(seq, file_name='out.pdf', snips=snips, scores=scores, fold=fold)


if __name__ == '__main__':
    args = docopt(__doc__, version='RaSE v1.0')
    main(args)
