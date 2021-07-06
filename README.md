# PopCanyon
**PopCanyon** discovers the under-methylated regions (UMRs) by solving Beta-binomial HMM on raw counts from multiple methylation files. In addition, based on previous detected UMRs, it can identify the differential methylated regions (DMRs) using Beta-binomial likelihood ratio test.
## Authors
- Xiaodong Cui (Xiaodong.Cui@bcm.edu)
- Jiejun Shi (jiejuns@uci.edu)
- Jianfeng Xu (jianfenx@bcm.edu)
- Wei Li (wei.li@uci.edu)
## Dependencies
- Python3 with following packages
  - numpy
  - scipy
  - pandas
  - copy
  - rpy2 (optional)
## Installation
No installation needed.
## Usage
There are two executable scripts: `./scripts/popcanyon.py` and `./scripts/MergeUMR.py`. `./scripts/utils.py` contains the functions required by `./scripts/popcanyon.py`.

`./scripts/popcanyon.py` includes the main functions for UMR and DMR calling. A detailed sub-commands help is available by typing:

```
$ python ./scripts/popcanyon.py UMR -h
usage: popcanyon UMR file1.txt file2.txt

Required arguments:
  --file FILE1 FILE2 FILE3 [FILE1 FILE2 FILE3 ...], -F FILE1 FILE2 FILE3 [FILE1 FILE2 FILE3 ...]
                        List of methylation files separated by spaces.
                        (default: None)

Optional arguments:
  --outdir OUTDIR, -O OUTDIR
                        Output directory for the results (default: ./result/)
  --use_cols USE_COLS [USE_COLS ...]
                        The index of the columns which should denote chr, pos,
                        strand, methy_reads, total_reads. Note: columns start
                        with 0 index (default: [0, 1, 2, 6, 7])
  --suffix SUFFIX       The suffix for the output filename (default: )
  --window_len WINDOW_LEN
                        Number of methylation sites for sliding (default: 3)
  --pval PVAL           p value for detecting the umrs (default: 0.01)
  --umr_len UMR_LEN     Minimum number of cpg sites in an umr (default: 10)
  --umr_fdr UMR_FDR     threshold for fdr value for detecting the umrs
                        (default: 0.05)
  --min_distance MIN_DISTANCE
                        Detected UMR length should overpass this threshold(Not
                        Used) (default: 2000.0)
  --min_sites MIN_SITES
                        Total number of methylation sites included in the UMR
                        region should overpass this threshold (default: 3)
  --wig                 Generate wig files for visualization (default: False)
  --viterbi_disable     Disable viterbi decoding process (default: False)
  --by_each             Besides popcanyon reporting UMR with all samples,
                        popcanyon will additionally report the umr according
                        to each sample. (default: False)
  --Nit NIT             The number of the iterations for Beta-Binomial HMM to
                        run (default: 6)
  --quiet               Print the processing messages (default: False)
  -h, --help            show this help and exit
```
	$ python ./scripts/popcanyon.py DMR -h

`./scripts/MergeUMR.py` is a downstream tool of UMR calling. It is used to merge the adjacent UMR if they meet specified threshold. A detailed sub-commands help is available by typing:

	$ python ./scripts/MergeUMR.py MergeUMR -h
