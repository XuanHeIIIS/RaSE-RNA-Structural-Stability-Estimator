# RaSE
RaSE - RNA structurAl Stability Estimator

<p align="center"><img src="img/trna.png"></p>

## Example

```
echo 'AACAGUCGAUCGAUCGAC' | ./code/RaSE.py --draw | tee out.txt && open out.pdf
```

## Help

```
RaSE - RNA structurAl Stability Estimator.

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
```
