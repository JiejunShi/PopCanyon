# PopCanyon
**PopCanyon** discovers the under-methylated regions (UMRs) by solving Beta-binomial HMM on raw counts from multiple methylation files. In addition, it can identify the differential methylated regions (DMRs) using Beta-binomial likelihood ratio test.
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

```
$ python ./scripts/popcanyon.py DMR -h
usage: popcanyon.py DMR --bed selection.bed --file file1.txt file2.txt -L treat ctrl

Required arguments:
  --file FILE1 FILE2 FILE3 [FILE1 FILE2 FILE3 ...], -F FILE1 FILE2 FILE3 [FILE1 FILE2 FILE3 ...]
                        List of methylation files separated by spaces.
                        (default: None)
  --labels treat treat ctrl [treat treat ctrl ...], -L treat treat ctrl [treat treat ctrl ...]
                        List of labels to denote the condition of the
                        methylation files; please input ``treat`` or ``ctrl``
                        and label numbers should be exactly the same as input
                        files (default: None)
  --bed FILE1.bed FILE2.bed [FILE1.bed FILE2.bed ...]
                        Do the DMR analysis only on the regions specified in
                        the bed file. (default: None)

Optional arguments:
  --outdir OUTDIR, -O OUTDIR
                        Output directory for the results (default: ./result/)
  --use_cols USE_COLS [USE_COLS ...]
                        The index of the columns which should denote chr, pos,
                        strand, methy_reads, total_reads. Note: columns start
                        with 0 index (default: [0, 1, 2, 6, 7])
  --suffix SUFFIX       The suffix for the output filename (default: )
  --dmr_fdr DMR_FDR     The threshold of fdr for reporting the dmr (default:
                        0.05)
  --quiet               Print the processing messages (default: False)
  -h, --help            show this help and exit
```

`./scripts/MergeUMR.py` is a downstream tool of UMR calling. It is used to merge the adjacent UMR if they meet specified threshold. A detailed sub-commands help is available by typing:

```
$ python ./scripts/MergeUMR.py MergeUMR -h
usage: python MergeUMR.py MergeUMR <positional arguments> [optional arguments]

Merge Adjacent UMRs if they meet specified threshold.

positional arguments:

  UMR                   UMRs called by PopCanyon
  MethRatio             MethRatio file from BSMAP which was used to call the
                        UMR

optional arguments:
  -h, --help            show this help message and exit
  -r , --methratio_max
                        Max methratio for UMR detection and Adjacent UMR
                        merging (default: 0.1)
  -a , --adjacent_max   Max distance of merged adjacent UMR (default: 500)
  -p , --prefix         Prefix of output file (default: out)
```

