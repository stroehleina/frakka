# frakka
Filter Kraken output files and calculate read-level and summary confidence score metrics per classified species.

Setting an appropriate Kraken2 [confidence score](https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring) threshold for species-level classification may appear rather arbitrary and will affect detection sensitivity. To take the guess work out of this, `frakka` creates per-species confidence 
score distribution metrics and plots.

# Installation

## GitHub

`git clone https://github.com/stroehleina/frakka`

## Conda

`TODO conda create -n frakka -c bioconda frakka`

# Usage

`python frakka.py --help`

```
usage: frakka.py [-h] [--version] [--kreport KREPORT] [--kout KOUT] [--species SPECIES] [--fof FOF] [--score SCORE] [--counts] [--taxid] [--plot] [--directory DIRECTORY] [--prefix PREFIX]
                 [--delim DELIM]

frakka - a tool to filter Kraken output files and calculate read-level and summary confidence score metrics per classified species

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit
  --kreport KREPORT, -k KREPORT
                        Kraken2 report file(s) (comma-separated if multiple, required) (default: None)
  --kout KOUT, -o KOUT  Kraken2 output file(s) (comma-separated if multiple, required, must be in same order as --kreport / -k) (default: None)
  --species SPECIES, -sp SPECIES
                        List of corresponding true / known species names or taxids (comma-separated if multiple, optional, must be in same order as -k and -o files, single species provided assumes
                        all files are same true / known species) (default: None)
  --fof FOF, -f FOF     Tab-separated file of one (-k, -o, -sp)-tuple per line (default: None)
  --score SCORE, -s SCORE
                        Confidence score threshold, only reads / counts higher than this score are reported (default: 0)
  --counts, -c          Report total counts per species with score > --score / -s instead of per-read reporting (default: False)
  --taxid, -t           Species input and output are NCBI taxIDs instead of species names (default: False)
  --plot, -p            NOT IMPLEMENTED Plot distribution of score per species (default: False)
  --directory DIRECTORY, -d DIRECTORY
                        NOT IMPLEMENTED Specify output directory (default: .)
  --prefix PREFIX, -x PREFIX
                        NOT IMPLEMENTED Specify output file prefix (default: .)
  --delim DELIM, -del DELIM
                        Specify output file delimiter (default: )
```

# Output

Tab-separated output to `STDOUT`:

`file` | `read_id` | `true_spec_taxid` | `K_spec_taxid` | `score`
--- | --- | --- | --- | --- 
file_1 | read_1 | species_A | species_A | 0.6
file_1 | read_2 | species_A | species_B | 0.4
file_1 | read_3 | species_A | species_A | 0.6
file_2 | read_a | species_C | species_C | 0.6
... | ... | ... | ... 
file_n | read_m | species_X | species_Y | 0.01

Tab-separated output to `STDOUT`, when `--counts` (and optionally, `--score`) are specified:

`file` | `K_spec_taxid` | `read_count`<br />(of reads above `--score`) | `median_score`<br />(of reads above `--score`)
--- | --- | --- | --- 
file_1 | species_A | count_1A | 0.6
file_1 | species_B | count_1B | 0.4
file_1 | species_C | count_1C | 0.6
file_2 | species_A | count_2A | 0.6
... | ... | ... | ... 
file_n | species_X | count_nX | 0.01

Logfile goes to `STDERR`.

## Graphical output

`TODO`

Plot distribution of score per species.

## Shiny confidence score sweep

`TODO`

Shiny app to allow to adjust confidence score and see how this affects number of reads binned into (selectable set of) species.

# Etymology

Frakka is partly "filter Kraken", partly [Old Norse for "spear"](https://en.wiktionary.org/wiki/frakkar) and partly inspired by the Torstyverse suite of tools by [Torsten Seemann](https://github.com/tseemann) which follow a similar naming pattern.
