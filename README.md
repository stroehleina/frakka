# Frakka
Filter Kraken output files and calculate read-level and summary confidence score metrics per classified species.

Setting an appropriate Kraken2 [confidence score](https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring) threshold for species-level classification may appear rather arbitrary and will affect detection sensitivity. To take the guess work out of this, Frakka creates per-species confidence score distribution metrics and plots.

`frakka --kreport [k2.tab1,k2.tab2,k2.tab3,...] --kout [k2.out1,k2.out2,k2.out3,...] [OPTIONS]`

```
--kreport (-k) [] - List of Kraken2 report files (required)
--kout (-o) [] - List of corresponding Kraken2 output files (required)
--fof (-f) - tab-separated file of report and output file pairs
--counts (-c) - report total counts for each species instead of per-read output
--score (-s) - confidence score threshold [0]
--species (-sp) [] - List of true species taxids / names
```
