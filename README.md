# Frakka
Filter Kraken output files and calculate read-level and summary confidence score metrics per classified species.

Setting an appropriate Kraken2 [confidence score](https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring) threshold for species-level classification may appear rather arbitrary and will affect detection sensitivity. To take the guess work out of this, Frakka creates per-species confidence 
score distribution metrics and plots.

# Installation

`TODO`

# Usage

`frakka --kreport [k2.tab1,k2.tab2,k2.tab3,...] --kout [k2.out1,k2.out2,k2.out3,...] [OPTIONS]`

```
--kreport (-k) [] - Comma-separated list of Kraken2 report files (required)
--kout (-o) [] -  Comma-separated list of corresponding Kraken2 output files (required)
--species (-sp) [] -  Comma-separated list of true species taxids / names (optional)
--fof (-f) FILE - tab-separated file of report and output file pairs (and optionally, a third column with species identifiers)
--counts (-c) BOOL - report total counts per species with score greater than the specified threshold (--score/-s) instead of per-read output
--score (-s) FLOAT - confidence score threshold [0] (optional)
--taxid (-t) BOOL - Species input and output are NCBI taxID instead of Species name
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

`file` | `K_spec_taxid` | `read_count` (of reads above `--score`) | `median_score` (of reads above `--score`)
--- | --- | --- | --- 
file_1 | species_A | count_1A | 0.6
file_1 | species_B | count_1B | 0.4
file_1 | species_C | count_1C | 0.6
file_2 | species_A | count_2A | 0.6
... | ... | ... | ... 
file_n | species_X | count_nX | 0.01

Logfile goes to `STDERR`.

# Graphical output

`TODO`

# Shiny confidence score sweep

`TODO`

# Etymology

Frakka is partly "filter Kraken", partly [Old Norse for "spear"](https://en.wiktionary.org/wiki/frakkar) and partly inspired by the Torstyverse suite of tools by [Torsten Seemann](https://github.com/tseemann) which follow a similar naming pattern.
