# PopCanyon
**PopCanyon** discovers the under-methylated regions (UMRs) by solving Beta-binomial HMM on raw counts from multiple methylation files. In addition, based on previous detected UMRs, it can identify the differential methylated regions (DMRs) using Beta-binomial likelihood ratio test.
## Authors
- Xiaodong Cui (Xiaodong.Cui@bcm.edu)
- Jiejun Shi (jiejuns@uci.edu)
- Wei Li (wei.li@uci.edu)
## Dependencies
- Python3 with following packages
  - numpy
  - pandas
  - copy
## Installation
No installation needed.
## Usage
There are two executable scripts: `./scripts/popcanyon.py` and `./scripts/MergeUMR.py`. `./scripts/utils.py` contains the functions required by `./scripts/popcanyon.py`.

`./scripts/popcanyon.py` includes the main functions for UMR and DMR calling. A detailed sub-commands help is available by typing:

	$ python ./scripts/popcanyon.py UMR -h
	$ python ./scripts/popcanyon.py DMR -h

`./scripts/MergeUMR.py` is a downstream tool for UMR calling. It is used to merge the adjacent UMR if they meet specified threshold. A detailed sub-commands help is available by typing:

	$ python ./scripts/MergeUMR.py MergeUMR -h
