#!/usr/bin/env python


"""RaSE - RNA structurAl Stability Estimate.

Compute stability.

Version: 1.1
Author: Fabrizio Costa [costa@informatik.uni-freiburg.de]

Usage:
  RaSE [-i <sequence>]
       [-k N] [-c N, --complexity=N] [-n N, --nbits=N] [-w N, --window_size=N]
       [-b N, --max_bp_span=N] [-p N, --avg_bp_prob_cutoff=N]
       [-r N, --hard_threshold=N] [-e N, --max_num_edges=N]
       [-l, --no_lonely_bps] [-t, --no_nesting]
       [--draw] [--jpg | --svg | --png | --pdf]
       [--verbose]
  RaSE (-h | --help | --version)

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
  --jpg                             Save images in jpg format.
  --svg                             Save images in svg format.
  --png                             Save images in png format.
  --pdf                             Save images in pdf format.
  -h --help                         Show this screen.
  --version                         Show version.
  --verbose                         Print more text.


"""
import sys
import numpy as np
import math
from docopt import docopt
from toolz import curry, compose
from toolz.curried import map
import matplotlib.pyplot as plt

from sklearn.metrics.pairwise import pairwise_kernels

from eden.graph import Vectorizer
from eden.display import draw_graph
from eden.display import draw_graph_set
from eden_rna.rnaplfold import fold as fold_rnaplfold
from eden_rna.rnafold import rnafold_wrapper as fold_rnafold
from eden.util import configure_logging, serialize_dict

import logging
logger = logging.getLogger(__name__)


def _save(text, full_out_file_name):
    with open(full_out_file_name, 'w') as f:
        for line in text:
            f.write("%s\n" % line.encode('utf8').strip())


def dotbracket(seq):
    """Compute the MFE in dotbracket notation."""
    seq_info, seq_struct = fold_rnafold(seq)
    return seq_struct


def _replace(seq, position, character):
    tokens = list(seq)
    tokens[position] = character
    return ''.join(tokens)


def compute_mfes(seq, mutations):
    """Given a sequence and the alternative sequence yield all mfes."""
    for position, character in enumerate(mutations):
        alt_seq = _replace(seq, position, character)
        yield dotbracket(alt_seq)


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
    fold = curry(fold_rnaplfold)(window_size=window_size,
                                 max_bp_span=max_bp_span,
                                 avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                                 hard_threshold=hard_threshold,
                                 max_num_edges=max_num_edges,
                                 no_lonely_bps=no_lonely_bps,
                                 nesting=nesting)
    return fold


def _window_reweight(boundaries, original_graph):
    graph = original_graph.copy()
    if boundaries is not None:
        begin, end = boundaries
        for i, u in enumerate(graph.nodes()):
            if begin <= i <= end:
                w = 1
            else:
                w = 0
            graph.node[u]['weight'] = w
    else:
        pass
    return graph


def make_fold_vectorize(complexity=3,
                        nbits=15,
                        fold=None,
                        boundaries=None):
    """Curry parameters in vectorizer."""
    vec = Vectorizer(complexity=complexity, nbits=nbits)
    vectorize = curry(lambda vec, graphs: vec.transform(graphs))(vec)

    cwindow_reweight = curry(_window_reweight)(boundaries)
    fold_vectorize = compose(vectorize, map(cwindow_reweight), fold)
    return fold_vectorize


def _make_string(gen):
    return ''.join(list(gen))


def _make_variation(seq, index, mutation):
    alternative = list(seq)
    alternative[index] = mutation
    return _make_string(alternative)


def _make_variations(seq, index=0, alphabet='ACGU'):
    for nt in alphabet:
        alternative = _make_variation(seq, index, nt)
        yield (nt, alternative)


def compute_stability(seq, alphabet='ACGU', fold_vectorize=None):
    """Compute the structural effects of single nt change.

    Specifically, compute the least similarity of the structure obtained
    by replacing each nucleotide with all possible mutations.
    """
    # TODO: parallelize indices
    matrix = np.zeros((len(alphabet), len(seq)))
    scores = list()
    mutations = list()
    vec = fold_vectorize([('header_placeholder', seq)])
    for index in range(len(seq)):
        _cmake_variations = curry(_make_variations)(index=index,
                                                    alphabet=alphabet)
        make_vecs = compose(fold_vectorize, _cmake_variations)
        variation_vecs = make_vecs(seq)
        sims = pairwise_kernels(vec, variation_vecs, metric='cosine')
        sims = sims.reshape((len(alphabet), 1))
        sims_variations = zip(sims, _cmake_variations(seq))
        score, alternative = min([(sim[0], variation[0])
                                  for sim, variation in sims_variations])
        scores.append(score)
        mutations.append(alternative)
        for i in range(4):
            matrix[i, index] = sims[i]
    return scores, mutations, matrix


def stability(seq, alphabet='ACGU', fold_vectorize=None):
    """Compute the structural effects of single nt change.

    Specifically, compute the least similarity of the structure obtained
    by replacing each nucleotide with all possible mutations.
    This function wraps compute_stability and post processes its output.
    """
    scores, mutations, matrix = compute_stability(
        seq, alphabet, fold_vectorize)
    mutations = _make_string(mutations)
    return mutations, scores, matrix


def stability_score(seq, pos=None, mutation=None, fold_vectorize=None):
    """stability."""
    vec1 = fold_vectorize([('header_placeholder', seq)])
    alternative = _make_variation(seq, pos, mutation)
    vec2 = fold_vectorize([('header_placeholder', alternative)])
    norm1 = vec1.dot(vec1.T)
    norm2 = vec2.dot(vec2.T)
    score = vec1.dot(vec2.T) / np.sqrt(norm1 * norm2)
    score = score[0, 0]
    return score


def serialize(seq, mutations, scores, k=5):
    """Pretty print of the alternative information and the relative scores."""
    num_digits = int(math.log10(len(seq))) + 1
    alt_tuples = zip(mutations, scores, seq)
    tuples = sorted(
        [(score, i, nt, alternative)
         for i, (alternative, score, nt) in enumerate(alt_tuples)])
    score_th = tuples[k][0]
    mfes = compute_mfes(seq, mutations)
    tuples = zip(mutations, scores, seq, mfes)
    spacer = '            '
    seq_str = spacer + '%s' % seq
    yield seq_str
    struct_str = spacer + '%s' % dotbracket(seq)
    yield struct_str
    for i, (alternative, score, nt, struct) in enumerate(tuples):
        if score < score_th:
            mark = '*'
        else:
            mark = ''
        if num_digits < 3:
            line = '%s %02d %s ' % (nt, i, alternative)
        elif num_digits == 3:
            line = '%s %03d %s ' % (nt, i, alternative)
        else:
            line = '%s %d %s ' % (nt, i, alternative)
        yield line + '%.2f %s %s' % (score, struct, mark)


class StructuralStabilityEstimator(object):
    """Class for the computation of the structural effects of nt change."""

    def __init__(self,
                 seq,
                 importance_window_begin=None,
                 importance_window_end=None,
                 alphabet='ACGU',
                 k=5,
                 complexity=3,
                 nbits=15,
                 window_size=150,
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
        self.seq = seq
        self.importance_window_begin = importance_window_begin
        self.importance_window_end = importance_window_end
        self.k = k
        self.complexity = complexity
        self.nbits = nbits
        self.mutations = None
        self.scores = None
        self.matrix = None
        window_size = max(window_size, len(seq))
        max_bp_span = max(max_bp_span, len(seq))
        self.fold = make_fold(window_size=window_size,
                              max_bp_span=max_bp_span,
                              avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                              hard_threshold=hard_threshold,
                              max_num_edges=max_num_edges,
                              no_lonely_bps=no_lonely_bps,
                              nesting=nesting)
        if self.importance_window_begin is not None \
                and self.importance_window_end is not None:
            begin = max(0, self.importance_window_begin)
            end = min(len(self.seq), self.importance_window_end)
            boundaries = (begin, end)
        else:
            boundaries = None
        self.fold_vectorize = make_fold_vectorize(complexity=self.complexity,
                                                  nbits=self.nbits,
                                                  fold=self.fold,
                                                  boundaries=boundaries)

    def stability_score(self, pos=None, mutation=None,
                        importance_semi_window=None):
        """stability."""
        if importance_semi_window is not None:
            begin = max(0, pos - importance_semi_window)
            end = min(len(self.seq), pos + importance_semi_window)
            boundaries = (begin, end)
        else:
            if self.importance_window_begin is not None \
                    and self.importance_window_end is not None:
                begin = max(0, self.importance_window_begin)
                end = min(len(self.seq), self.importance_window_end)
                boundaries = (begin, end)
            else:
                boundaries = None
        self.fold_vectorize = make_fold_vectorize(complexity=self.complexity,
                                                  nbits=self.nbits,
                                                  fold=self.fold,
                                                  boundaries=boundaries)
        return stability_score(self.seq,
                               pos=pos,
                               mutation=mutation,
                               fold_vectorize=self.fold_vectorize)

    def compute_stability_scores(self):
        """Computation of the structural effects of nt change."""
        self.mutations, self.scores, self.matrix = stability(
            self.seq,
            alphabet='ACGU',
            fold_vectorize=self.fold_vectorize)

    def output(self):
        """Pretty print the computed results."""
        return serialize(self.seq,
                         self.mutations,
                         self.scores,
                         k=self.k)

    def _update_graph(self, graph, scores, mutations):
        for u in graph.nodes():
            pos = graph.node[u]['position']
            graph.node[u]['stability'] = scores[pos]
            graph.node[u]['alternative'] = mutations[pos]
        return graph

    def _draw(self,
              seq,
              file_name=None,
              mutations=None,
              scores=None,
              fold=None):
        seqs = [('header_placeholder', seq)]
        graphs = fold(seqs)
        graph = graphs.next()
        graph = self._update_graph(graph, scores, mutations)

        opts = {'size': 14, 'font_size': 9, 'layout': 'KK',
                'vertex_border': False, 'vertex_size': 400,
                'vertex_alpha': 0.5, 'invert_colormap': True,
                'secondary_vertex_label': 'alternative',
                'vertex_color': 'stability', 'colormap': 'OrRd',
                'edge_label': None, 'edge_alpha': 0.2,
                'dark_edge_alpha': 0.8,
                'ignore_for_layout': 'nesting'}
        draw_graph(graph, file_name=file_name, **opts)

    def draw(self, file_name=None):
        """Graphical representation of the RNA folded structure.

        Nodes encode the original nt and the most de-stabilizing alternative.
        """
        self._draw(self.seq,
                   file_name=file_name,
                   mutations=self.mutations,
                   scores=self.scores,
                   fold=self.fold)

    def plot_instability(self, file_name=None):
        """Graph of unstability in the RNA sequence."""
        size = len(self.scores) / 2.5
        fig = plt.figure(figsize=(size, 2.5))
        ax1 = fig.add_subplot(111)
        width = 1
        x = np.array(range(len(self.scores))) - width / 2.0
        y = 1 - np.array(self.scores)
        ax1.bar(x, y, width, alpha=0.3)
        ax1.set_xticks(range(len(self.seq)))
        ax1.set_xticklabels(list(self.seq))
        ax1.set_xlim(-1, len(self.seq))
        ax1.set_ylim(0, 1)
        ax2 = ax1.twiny()
        ax2.set_xlim(-1, len(self.seq))
        ax2.set_xticks(range(len(self.seq)))
        ax3 = ax1.twiny()
        ax3.yaxis.set_visible(False)
        ax3.set_xlim(-1, len(self.seq))
        ax3.set_xticks(range(len(self.seq)))
        ax3.set_xticklabels(list(self.mutations))
        pos = ax3.get_position()
        pos = [pos.x0, pos.y0, pos.width, 0.001]
        ax3.set_position(pos)
        if file_name is None:
            plt.show()
        else:
            plt.savefig(file_name,
                        bbox_inches='tight',
                        transparent=True,
                        pad_inches=0)
            plt.close()

    def plot_stability_matrix(self, file_name=None):
        """Plot matrix of unstability in the RNA sequence."""
        size = len(self.seq) / 2.5
        plt.figure(figsize=(size, 2.5))
        plt.imshow(self.matrix,
                   interpolation='none',
                   cmap=plt.get_cmap('YlOrRd'))
        plt.yticks(range(4), ['A', 'C', 'G', 'U'], fontsize=12)
        plt.xticks(range(len(self.seq)), fontsize=12)
        if file_name is None:
            plt.show()
        else:
            plt.savefig(file_name,
                        bbox_inches='tight',
                        transparent=True,
                        pad_inches=0)
            plt.close()

    def draw_all(self, file_name=None):
        """Plot all k most unstable structures."""
        infos = zip(self.mutations, self.scores, self.seq)
        tuples = sorted(
            [(score, i, nt, alternative)
             for i, (alternative, score, nt) in enumerate(infos)])
        tuples = tuples[:self.k]
        graphs = []
        for score, position, nt, alternative in tuples:
            header = '%s %d %s' % (nt, position, alternative)
            alt_seq = _replace(self.seq, position, alternative)
            _graphs = self.fold([(header, alt_seq)])
            graph = _graphs.next()
            graphs.append(graph)

        opts = {'size': 10, 'font_size': 9, 'colormap': 'rainbow',
                'vertex_border': False, 'vertex_size': 200,
                'vertex_alpha': 0.4, 'vertex_color': '_label_',
                'edge_alpha': 0.2, 'edge_label': None,
                'dark_edge_alpha': 0.8,
                'ignore_for_layout': 'nesting', 'layout': 'KK',
                'n_graphs_per_line': 3}
        draw_graph_set(graphs, file_name=file_name, **opts)


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
    draw = args['--draw']
    jpg = args['--jpg']
    svg = args['--svg']
    png = args['--png']
    pdf = args['--pdf']

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
    rase = StructuralStabilityEstimator(seq,
                                        alphabet='ACGU',
                                        k=k,
                                        complexity=complexity,
                                        nbits=nbits,
                                        window_size=window_size,
                                        max_bp_span=max_bp_span,
                                        avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                                        hard_threshold=hard_threshold,
                                        max_num_edges=max_num_edges,
                                        no_lonely_bps=no_lonely_bps,
                                        nesting=nesting)
    # print: nt pos, original nt, most de-stabilizing nt, dotbracket, score
    for line in rase.transform(seq):
        print(line)

    # if drawing is required use the folding algorithm to compute the graph
    if draw:
        suffix = 'pdf'
        if jpg:
            suffix = 'jpg'
        if svg:
            suffix = 'svg'
        if png:
            suffix = 'png'
        if pdf:
            suffix = 'pdf'
        structure_fname = 'structure.' + suffix
        score_fname = 'score.' + suffix
        all_plots_fname = 'structures.' + suffix
        rase.draw(file_name=structure_fname)
        rase.plot(file_name=score_fname)
        rase.draw_all(file_name=all_plots_fname)


if __name__ == '__main__':
    args = docopt(__doc__, version='RaSE v1.1')
    main(args)
