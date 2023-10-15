# Augmented Thresholds for MONI
See PHONI, which this method is a variant of: https://github.com/koeppl/phoni

Modified version of one-pass MONI using runs compressed BWT index to generate matching statistics for a pattern. Uses threshold LCE trade-off as described in [Augmented Thresholds for MONI](10.1109/DCC55655.2023.00035). Cite this paper if you use this tool

Requires the same dependencies as PHONI to run all scripts/benchmarks: we borrow their description below modified for differences in running this version.
## Preparations

We require the pattern and the text to be available in form of sequences stored in the `.fa` (FASTA) format.
To use our solution, you need to have recent `cmake`, `g++`, `zsh`, and `python 3` installed.

We need the following python 3 packages for extracting and concatenating `.fa` files:
```console
	pip3 install biopython
	pip3 install fastaparser
	pip3 install psutil
```


```console
git clone --branch phoni https://github.com/drnatebrown/aug_phoni
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

### Building the index

To *build* the index we use the command `aug build` from the build directory. Alternatively, replace `aug` with `phoni` to build the version without augmented thresholds.

``` console
aug build \
-r <filename of the reference> \
-t <number of threads> \
-g <grammar format> \
-f <input file is a fasta file> \
```
For example, to build the aug-phoni index for the file `yeast.fa` using 4 `threads` and the `plain` grammar we run from the `build` folder:
``` conole
python3 aug build -r ../data/yeast.fa -f -t 4 -g plain
```

This command will produce `yeast.fa.aug` and `yeast.fa.plain.slp` in the `data` folder, which represent the `aug-phoni` index.

### Querying the index

To *query* the index we use the command `aug query` from the build directory. We replace `aug` with `phoni` to instead run without augmented, assuming PHONI was also built.

``` console
aug ms \
-i <filename of the reference> \
-p <fasta pattern file> \
-g <grammar format> \
```
For example, to query the phoni index for the file `yeast.fa` using the `plain` grammar with the pattern `samples.fa` we run from the `build` folder:
``` console
python3 aug ms -r ../data/yeast.fa -p ../data/samples.fa -g plain
```

This command will produce `samples.fa.positions` and `samples.fa.lengths` in the `data` folder, which represent the matching staistics *positions* and *lengths* of `samples.fa` against `yeast.fa`, respectively.

### Benchmarks

We provide a script and benchmark files to evaluate Aug PHONI described in the paper:

C. Martínez-Guardiola, N. K. Brown, F. Silva-Coira, D. Köppl, T. Gagie and S. Ladra, "Augmented Thresholds for MONI," 2023 Data Compression Conference (DCC), Snowbird, UT, USA, 2023, pp. 268-277, doi: 10.1109/DCC55655.2023.00035.

In our experiments we used the file
 - [chr19.1000.fa.xz](http://dolomit.cs.tu-dortmund.de/tudocomp/chr19.1000.fa.xz) as our text dataset, and prefixes of it as our pattern

We have a shell script `benchmark.sh` for an automatic benchmark.
For this to work, some variables in it has to be set, meaning it is necessary to download and compile those projects individually, and the set the corresponding variables in `benchmark.sh` manually. Finally, the output of `benchmark.sh` can be processed by [sqlplots](https://github.com/koeppl/sqlplot) to generate the plots shown in the paper.
