[![DOI](https://zenodo.org/badge/74471783.svg)](https://zenodo.org/badge/latestdoi/74471783)

# RaSE - RNA structurAl Stability Estimator
RaSE is a python program that uses the graph vectorization technique of [EDeN](https://github.com/fabriziocosta/EDeN) to compute a score indicative of the structural stability responsibility of each nucleotide in an RNA sequence. The score is computed as the similarity of the structure obtained by changing a single nucleotide with respect to the original structure. Out of the 3 possible alternatives, only the one which yields the largest difference is reported. The problem of assessing the similarity between two structures is cast in the equivalent problem of assessing the similarity between two graphs which encode the structural information. Structures are computed using the [RNAplfold](https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html) program. The graph embedding represents nucleotides as nodes labeled with the one letter code (A|C|G|U); relations between backbone bounds and base pairs with a probability higher than --hard_threshold are repersented as conjunctive edges, relations between base pairs with probability in the interval --avg_bp_prob_cutoff and --hard_threshold are represented as disjunctive edges. Graphs are transformed into a high dimensional saprse vector representation using the graph vectorization technique of [EDeN](https://github.com/fabriziocosta/EDeN). The similarity between the corresponding vectors is then computed as the [cosine similarity](https://en.wikipedia.org/wiki/Cosine_similarity).


## Input
RaSe takes in input a RNA sequence as a single string of one letter code (A|C|G|U). The string can be provided via the stdin or via the -i flag.


## Output
RaSe outputs to stdout a space separated tabular file with the following format: the first line contains the Minimum Free Energy structure (MFE) in dotbracket notation; the following lines contain the nucleotide (nt) position, the nt one letter code (A|C|G|U), the nt code for the mutation that most changes the computed structure, the similarity score between the original structure and the structure obtained by the mutation, the MFE of the mutated sequence and an optional character * marking the top dissimilar cases.  

```
             (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))).
  0 G C 0.63 ((((.((((((.((.......)).)))))).)))).......((((((...))))))(((((.......)))))........ 
  1 C G 0.23 (((((((.((((((((.........((((((.......))))))..))))))))...(((((.......)))))))))))). *
```
RaSE can optionally produce image files in various formats (jpg, png, svg, pdf). When invoked with the flag --draw the following files are produced: structure.[format], plot.[format], structures.[format] 

The structure image depicts the graph encoding of the most probable RNA structure: edges between backbone bounds and base pairs with a probability higher than --hard_threshold are displayed with a solid line, edges between base pairs with probability in the interval --avg_bp_prob_cutoff and --hard_threshold are displayed with a dashed line. The node label is composed of the original nt (above) and the mutation that most changes the computed structure (below). The color intensity is proportional to 1 - similarity, so that darker nodes are the ones that have the largest effect on the structure.

<p align="center"><img src="img/trna.png"></p>

The plot image depicts the nt position on the top x axis, the original nt on the bottom x axis, the mutation that most changes the computed structure on the bottom x axis but inside the plot, the score = 1 - similarity on the y axis, so that the highest bar corresponds to the mutation that has the largest effect on the structure.

<p align="center"><img src="img/plot.png"></p>

The structures image depicts the graph encoding of the individual k mutations that most changes the computed structure. The title associated with each graph is composed of the original nt, the position and the mutation that most changes the computed structure.

<p align="center"><img src="img/structures.png"></p>

## Library
RaSE exposes several functions that can be used inside other projects. See examples of use in the 
 [Jupyter notebook](https://github.com/fabriziocosta/RaSE/blob/master/RNAStructuralStabilityEstimator.ipynb)

## Example

```
echo 'GUGGACGUGCCGGAGUGGUUAUCGGGCAUGACUAGAAAUCAUGUGGGCUUUGCCCGCGCAGGUUCGAAUCCUGCCGUUCACG' | ./RaSE.py
             (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))).
  0 G C 0.63 ((((.((((((.((.......)).)))))).)))).......((((((...))))))(((((.......)))))........ 
  1 U C 0.61 ((((.((((((.((.......)).))))))...........(((((((...)))))))((((.......))))))))..... 
  2 G A 0.64 .(((.((((((.((.......)).)))))).)))........((((((...))))))(((((.......)))))........ 
  3 G U 0.70 (((..((((((.((.......)).))))))............((((((...))))))(((((.......)))))....))). 
  4 A G 0.35 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
  5 C A 0.64 (((((.(((((.((.......)).))))).............((((((...))))))(((((.......)))))..))))). 
  6 G A 0.69 (((((((((((.((.......)).))))))............((((((...))))))(((((.......)))))..))))). 
  7 U G 0.71 (((((((((((.((.......)).)))((((.......))))((((((...)))))).((((.......)))))))))))). 
  8 G A 0.48 (((((((..((.((.......)).))(((((.......)))))(((((...))))).(((((.......)))))))))))). 
  9 C G 0.23 (((((((.((((((((.........((((((.......))))))..))))))))...(((((.......)))))))))))). *
 10 C U 0.48 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 11 G A 0.55 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 12 G C 0.32 (((((((....((((.(((..((((((.((...........(((((((...))))))))).)))))))))))))))))))). *
 13 A C 0.37 (((((((.(((.((((((((.((((......)).)))))))))).))).........(((((.......)))))))))))). 
 14 G C 0.46 (((((((.((.(((((((((((.....)))))))).((((.(((((((...)))))))..))))....))).))))))))). 
 15 U G 0.95 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 16 G A 0.96 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 17 G A 0.96 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 18 U G 0.53 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 19 U C 0.47 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 20 A G 0.32 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 21 U A 0.37 (((((((((((.(.........).))))..............((((((...))))))(((((.......)))))))))))). 
 22 C A 0.43 (((((((.((.(((.(((((((.....)))))))(((....(((((((...)))))))....)))...))).))))))))). 
 23 G A 0.47 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 24 G C 0.45 (((((((.((.(((.(((((((.....)))))))(((....(((((((...)))))))....)))...))).))))))))). 
 25 G U 0.42 (((((((.((.(((.(((((((.....)))))))(((....(((((((...)))))))....)))...))).))))))))). 
 26 C A 0.45 (((((((.((.(((.(((((((.....)))))))(((....(((((((...)))))))....)))...))).))))))))). 
 27 A G 0.48 (((((((.(((.((.......)).)))((((.......))))((((((...))))))(((((.......)))))))))))). 
 28 U G 0.45 (((((((..((((........))))((((((...(((....(((((((...)))))))....)))...))))))))))))). 
 29 G C 0.78 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 30 A G 0.36 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 31 C U 0.76 (((((((((((.((.......)).))))(((((...))))).((((((...))))))(((((.......)))))))))))). 
 32 U C 0.40 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 33 A C 0.45 (((((((.......((((((.(((((.....)))))))))))((((((...))))))(((((.......)))))))))))). 
 34 G C 0.33 (((((((.(((.....)))...((((((...(((((.....)))))....)))))).(((((.......)))))))))))). 
 35 A G 0.46 (((((((.((.(((.(((((((.....))))))).(((((.(((((((...)))))))..)))))...))).))))))))). 
 36 A G 0.88 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 37 A U 0.47 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 38 U C 0.47 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 39 C G 0.46 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 40 A C 0.44 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 41 U C 0.50 (((((((((((.((.......)).)))).............(((((((...)))))))((((.......)))).))))))). 
 42 G U 0.35 (((((((.(((.....)))...((((((...(((((.....)))))....)))))).(((((.......)))))))))))). 
 43 U A 0.67 (((((((((((.((.......)).))))................((((...))))..(((((.......)))))))))))). 
 44 G U 0.41 (((((((.(((.....)))...((((((...((((........))))...)))))).(((((.......)))))))))))). 
 45 G C 0.30 (((((((...(((.((((....((.((((((.......)))))).))..))))))).(((((.......)))))))))))). *
 46 G C 0.40 (((((((.(((.....)))...((((((.(.(((.(.....).))).)..)))))).(((((.......)))))))))))). 
 47 C U 0.50 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 48 U G 0.94 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 49 U A 0.65 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 50 U C 0.95 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 51 G C 0.92 (((((((((((.((.......)).))))..............(((((.....)))))(((((.......)))))))))))). 
 52 C G 0.31 (((((((.((.(((.(((((((.....)))))))..((((.(((((((....))).))))))))....))).))))))))). *
 53 C A 0.16 (((((((((((.((.......)).))))..............(((.((...)).)))(((((.......)))))))))))). *
 54 C A 0.35 (((((((((((.((.......)).))))..............((.(((...))).))(((((.......)))))))))))). 
 55 G C 0.27 (((((((..((((........))))((((((.......))))))((((...))))..(((((.......)))))))))))). *
 56 C A 0.34 (((((((.(((.....)))..(((((((...(((.(.....).)))....)))))))(((((.......)))))))))))). 
 57 G U 0.46 (((((((.((.(((.(((((((.....)))))))(((.....((((((...)))))).....)))...))).))))))))). 
 58 C U 0.49 (((((((.((.(((.(((((((.....)))))))..((((((((((((...)))))))).))))....))).))))))))). 
 59 A C 0.35 (((((((.((((((.((((...((((((...(((.(.....).)))....)))))).)))).))))......))))))))). 
 60 G C 0.37 (((((((.((.(((.(((((((.....)))))))(((....(((((((...)))))))....)))...))).))))))))). 
 61 G C 0.36 (((((((.((((((((.((...((((((...(((.(.....).)))....)))))).)).))))))......))))))))). 
 62 U C 0.49 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 63 U C 0.53 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 64 C A 0.96 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 65 G C 0.59 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 66 A C 0.60 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 67 A C 0.36 (((((((.((.(((((.(..(((..((((((.......))))((((((...)))))))).))).).))))).))))))))). 
 68 U A 0.91 (((((((((((.((.......)).))))..............((((((...))))))(((((((...)))))))))))))). 
 69 C U 0.49 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))). 
 70 C G 0.47 (((((((..((((........))))((((((...(((....(((((((...)))))))....)))...))))))))))))). 
 71 U C 0.50 (((((((.((.(((.(((((((.....)))))))(((....(((((((...)))))))....)))...))).))))))))). 
 72 G U 0.50 (((((((((((.((.......)).)))).............(((((((...))))))).(((.......)))..))))))). 
 73 C A 0.70 (((((((((((.((.......)).)))).............(((((((...)))))))((((.......)))).))))))). 
 74 C G 0.68 (((((((((((.((.......)).))))).............((((((...))))))(((((.......))))).)))))). 
 75 G U 0.64 (((((((((((.((.......)).))))))............((((((...))))))(((((.......)))))..))))). 
 76 U G 0.37 .(((.((((((.((.......)).)))))).)))........((((((...))))))(((((.......)))))........ 
 77 U A 0.66 .(((.((((((.((.......)).)))))).)))........((((((...))))))(((((.......)))))........ 
 78 C A 0.64 .(((.((((((.((.......)).)))))).)))........((((((...))))))(((((.......)))))........ 
 79 A C 0.64 .(((.((((((.((.......)).)))))).)))........((((((...))))))(((((.......)))))........ 
 80 C G 0.68 .(((.((((((.((.......)).)))))).)))........((((((...))))))(((((.......)))))........ 
 81 G U 0.96 (((((((((((.((.......)).))))..............((((((...))))))(((((.......)))))))))))).
```

## Help

```
RaSE - RNA structurAl Stability Estimate.

Compute stability.

Version: 1.0
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

```
