# K-Near Neighbour Analysis (NN) of ECFP or GOBBI RDKit fingerprints
====================================================================

This is a script to perform near neighbour analysis of compound sets.
Can be used to calculate intra- or inter- similarity for compounds.
Intra-similarity is automatically performed if no second file provided.
Can take multiple files using * wildcard.

Main dependencies: 'RDKit'

NN_analysis.py [options]
-------------------------------

Options:
--------

```-h, --help```            show this help message and exit
  

```--i1=FILE```             First inputfile

```--i2=FILE```             Second inputfile (optional)

```-o FILE, --out=FILE```   Outputfile (optional)

```--d1=DELIM1, --delim1=DELIM1```     Delimiter for first file (default: ' ')

```--d2=DELIM2, --delim2=DELIM2```     Delimiter for second file (default: ' ')

```--delimcol1=DELIMCOL1```  Delimiter column of smiles in first file (default:0)

```--delimcol2=DELIMCOL2```  Delimiter column of smiles in second file (default:0)

```-b BITS, --bits=BITS```   No. bits in morgan fingerprint (default:2048)

```-r RADIUS, --radius=RADIUS```       Radius of morgan fingerprint (default:2)

```-t TOPN, --topn=TOPN```   Top n neighbors (default:1)

```-n NCORES, --ncores=NCORES```       No. cores (default:1)

```--pf1```                  Parallel process list of first files, not lines per file

```--pf2```                  Parallel process list of second files, not lines per file

```--gobbi```                Use Gobbi 2D pharmacophore FPs

```--mw=MW```                Max. Mw filter