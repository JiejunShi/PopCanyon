from utils import *
import argparse
import pandas as pd
import time
import sys
import os


############ 1. Parser Functions ############
def parse_arguments(args=None):
  parser = \
    argparse.ArgumentParser(
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""

  ``popcanyon`` discovers the under-methylated regions (UMRs) by solving Beta-binomial HMM on raw counts from multiple methylation files. In addition, based on previous detected UMRs, it can identify the differential methylated regions (DMRs) using Beta-binomial likelihood ratio test.
  The UMR and DMR analysis will be always performed for the entire samples. However, you can specify ``--by_each`` option to let the software do the same analysis per sample.

  If you are only interested in discovering UMRs, use the ``UMR`` mode instead.
  If you are interested in finding the DMRs for the previously discovered UMRs or regions, use the ``DMR`` mode.
  If you are interested in UMRs and DMRs, use the ``UMR-DMR`` mode.

  A detailed sub-commands help is available by typing:

  popcanyon UMR -h

  popcanyon DMR -h


  """,
          epilog='example usages:\n'
                 'popcanyon UMR -F file1.txt file2.txt --outdir ./results \n\n'
                 'popcanyon DMR --bed selection.bed --file file1.txt file2.txt -L treat ctrl\n\n'
                 ' \n\n',
          conflict_handler='resolve')

  parser.add_argument('--version', action='version',
                      version='popcanyon {}'.format(0.1))
  subparsers = parser.add_subparsers(
      title="commands",
      dest='command',
      description='subcommands',
      help='subcommands',
      metavar='')
  # bins mode options
  subparsers.add_parser(
      'UMR',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      parents=[conditional_args(case='UMR')
               ],
      help="The user provides methylation files and popcanyon will call the UMR region",
      add_help=False,
      usage='popcanyon UMR file1.txt file2.txt\n')
  # UMR-DMR mode arguments
  subparsers.add_parser(
      'UMR-DMR',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      parents=[conditional_args(case='UMR-DMR')
               ],
      help="The user provides methylation files with their corresponding labels: "
           "[treat or ctrl]. Popcanyon will report the region of UMRs and compare "
           "these UMRs to between two conditional sets of different samples",
      usage='popcanyon -F file1.txt file2.txt file3.txt -L treat treat ctrl\n',
      add_help=False)
  # DMR mode
  subparsers.add_parser(
      'DMR',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      parents=[conditional_args(case='DMR')
               ],
      help="The user provides a BED file that contains all regions "
           "that should be considered for the differential analysis. A "
           "common use is to compare the methylation changes between two conditional"
           " sets of different samples for a set of UMR regions.",
      usage='%(prog)s --bed selection.bed --file file1.txt file2.txt -L treat ctrl\n',
      add_help=False)
  return parser


def conditional_args(case='UMR'):
  parser = argparse.ArgumentParser(add_help=False)
  required = parser.add_argument_group('Required arguments')
  # define the arguments
  required.add_argument('--file', '-F',
                      metavar='FILE1 FILE2 FILE3',
                      help='List of methylation files separated by spaces.',
                      nargs='+',
                      required=True)
  optional = parser.add_argument_group('Optional arguments')

  optional.add_argument('--outdir', '-O', default='./result/',
  				help='Output directory for the results')

  optional.add_argument('--use_cols', default=[0,1,2,6,7], type=int, nargs='+',
          help='The index of the columns which should denote chr, pos, strand, methy_reads, total_reads. Note: columns start with 0 index')

  optional.add_argument('--suffix', default='', type=str, help='The suffix for the output filename')

  if case == 'UMR' or case == 'UMR-DMR':
    optional.add_argument('--window_len', default=3, type=int,
                          help='Number of methylation sites for sliding')

    optional.add_argument('--pval', default=0.01, type=float,
                          help='p value for detecting the umrs')

    optional.add_argument('--umr_len', default=10, type=int, help='Minimum number of cpg sites in an umr')

    optional.add_argument('--umr_fdr', default=0.05, type=float,
                          help='threshold for fdr value for detecting the umrs')

    # umr_parser.add_argument('--P_REGION_THRESHOLD', default=0.01, type=float,
    #                         help='p value for combining the sites together')
    optional.add_argument('--min_distance', default=2e3, type=int,
                          help='Detected UMR length should overpass this threshold(Not Used)')

    optional.add_argument('--min_sites', default=3, type=int,
                          help='Total number of methylation sites included in the UMR region should overpass this threshold')

    optional.add_argument('--wig', action='store_true',
                          help='Generate wig files for visualization')

    optional.add_argument('--viterbi_disable', action='store_true', help='Disable viterbi decoding process')

    optional.add_argument('--by_each', action='store_true',
     help='Besides popcanyon reporting UMR with all samples, '
      'popcanyon will additionally report the umr according to each sample.')

    optional.add_argument('--Nit', default=6, type=int, help='The number of the iterations for Beta-Binomial HMM to run')

	# parameters specified for UMR-DMR and DMR
  if case == 'DMR' or case == 'UMR-DMR':
    required.add_argument('--labels', '-L',
    	                  metavar='treat treat ctrl',
        	              help='List of labels to denote the condition of the methylation files; please input ``treat`` or ``ctrl`` and \
                        label numbers should be exactly the same as input files',
            	          nargs='+',
                	      required=True)
    optional.add_argument('--dmr_fdr',
                        help='The threshold of fdr for reporting the dmr',
                        default=0.05,
                        type=float)
	# parameters for DMR
  if case == 'DMR':
    required.add_argument('--bed',
                          help='Do the DMR analysis only on '
                          'the regions specified in the bed file.',
                          metavar='FILE1.bed FILE2.bed',
                          nargs='+',
                          required=True)

  # all other optional parameters
  optional.add_argument('--quiet', action='store_true',
                          help='Print the processing messages')
  optional.add_argument('-h', '--help', action='help', help='show this help and exit')
  return parser


# pre-check the inputs' requirements
def process_args(args=None):
  args = parse_arguments().parse_args(args)
  # add suffix to the final output filename
  if len(args.suffix) > 0:
    args.suffix = '.' + args.suffix

  # popcanyon will only read 5 columns of the input files
  if len(args.use_cols) != 5:
    raise IndexError('The total columns should be exactly 5!')

  # check the labels number and files number
  if args.command == 'DMR' or args.command == 'UMR-DMR':
    if len(args.labels) != len(args.file):
      raise NameError('The number of labels do not match number of files')

  # create the output folder
  check_and_create_dir(args.outdir)

  return args


############ 2. Main Functions ############
# print function
def disp(txt):
    if not args.quiet:
      if args.command == 'UMR':
        print('[UMR @{}] {}'.format(time.asctime(), txt), file=sys.stderr)
      if args.command == 'DMR':
        print('[DMR @{}] {}'.format(time.asctime(), txt), file=sys.stderr)
      if args.command == 'UMR-DMR':
        print('[UMR-DMR @{}] {}'.format(time.asctime(), txt), file=sys.stderr)

# main function
def main(args=None):
  # 2.1 preparation for reading the files
  sample_names = [file_name.split('/')[-1].split('.')[0] for file_name in args.file]
  id2sample = dict(zip(range(len(sample_names)), sample_names))
  chs = ['chr' + str(i) for i in range(1, 23)]  + ['chrX','chrY']
  # chs = ['chr1', 'chr2', 'chr3', 'chr7']
  n_files = len(sample_names) # number of files
  meth_files = read_methy_files(args.file, cols=args.use_cols) # read all the files into memory

  # 2.2 according to different mode and choose different modules
  # 2.2.1 if mode is not 'DMR', it will always call the umr first
  if args.command == 'UMR' or args.command == 'UMR-DMR':
    # if by_each is enabled, the UMR detection will first call UMR independently and then call UMRs for all
    if args.by_each and n_files > 1:
      umr_all = {i: {} for i in range(n_files)}
      umr_beds = {i: {} for i in range(n_files)}
      for ch in chs:
        # preparing for call all the meth profiles together
        # ch: chrome number
        # chr_sub: methy info for a particular chrome
        disp('Reading files for {} ...'.format(ch))
        # meth_files = read_methy_files(args.file, ch, cols=args.use_cols)
        for i in range(n_files):
          disp('Processing on {} in {}'.format(ch, id2sample[i]))
          meth_file = meth_files[i]
          chr_sub = meth_file[meth_file.chr == ch]
          # when chr_sub is not empty
          if chr_sub.shape[0] > 0:
            chr_sub = merge_strand(chr_sub)
            # downsample the reads proportionally on the sites with high duplicate reads
            chr_sub.loc[:,'methy'], chr_sub.loc[:,'unmethy'] = avoid_overflow_scaling(chr_sub.loc[:,'methy'].values, chr_sub.loc[:,'unmethy'].values)
            df_umr = call_umr_multi(chr_sub, ch, args)
            # filtering strategy
            umr_all[i][ch] = df_umr
            umr_beds[i][ch] = umr2bed(df_umr, which_fdr='umr_fdr')

        # write down all the info
        # modified: wrap up the combine each chromes into a big data frame code into a function 2017/12/29
      for i in range(n_files):
        disp('Writing results for {}'.format(id2sample[i]))
        # sample_all: one sample's umr detailed info
        # umr_all: includes all detailed info for umr
        df_chs_write_umr_all = combine_list_dfs_to_write(umr_all[i], chs)
        df_chs_write_umr_bed = combine_list_dfs_to_write(umr_beds[i], chs)
        # write the results
        df_chs_write_umr_all.to_csv('{}.umr.txt'.format(args.outdir+id2sample[i]+args.suffix),
                                    sep='\t', index=False)
        df_chs_write_umr_bed.to_csv('{}.umr.bed'.format(args.outdir+id2sample[i]+args.suffix),
                                    sep='\t', index=False,
                                    header=['# chr', 'chromStart', 'chromEnd', 'name', 'score',
                                            'strand', 'thickStart', 'thickEnd', 'itemRgb'])

    # 2.2.2 popcanyon will always report the UMR for all the samples combined
    #### UMR mode for all the samples
    umr_all = {}
    umr_beds = {}
    dmr_all = {}
    dmr_beds = {}
    for ch in chs:
      # preparing for call all the meth profiles together
      # ch: chrome number
      # chr_sub: methy info for a particular chrome
      # chr_subs: list of chr_sub
      disp('Processing on {}'.format(ch))
      chr_subs = []
      for i in range(n_files):
        meth_file = meth_files[i]
        chr_sub = meth_file[meth_file.chr == ch]
        chr_sub = merge_strand(chr_sub)
        # downsample the reads proportionally on the sites with high duplicate reads
        chr_sub.loc[:,'methy'], chr_sub.loc[:,'unmethy'] = avoid_overflow_scaling(chr_sub.loc[:,'methy'].values, chr_sub.loc[:,'unmethy'].values)
        chr_subs.append(chr_sub)
      # merge all the data in to a big data frame
      chr_merge = merge_methy_profiles(chr_subs)
      chr_merge = geometric_normalize(chr_merge)
      # when chr_sub is not empty
      if chr_merge.shape[0] > 0:
        df_umr_merge = call_umr_multi(chr_merge, ch, args)
        umr_all[ch] = df_umr_merge
        umr_beds[ch] = umr2bed(df_umr_merge, which_fdr='umr_fdr')

        #2.2.3 [SUB MODE] if UMR-DMR mode, then popcanyon will subsequently proceed dmr analysis
        if args.command == 'UMR-DMR':
          disp('Computing the DMR for {}'.format(ch))
          df_dmr_merge = call_dmr_multi(df_umr_merge, chr_merge, args)
          dmr_all[ch] = df_dmr_merge
          # umr2bed: convert the dataframe (umr or dmr) to bed format;
          # it is not just for umr
          dmr_beds[ch] = umr2bed(df_dmr_merge, which_fdr='dmr_fdr')


    # write the results for imr and umr
    disp('Writing the results')
    df_chs_write_umr_all = combine_list_dfs_to_write(umr_all, chs)
    df_chs_write_umr_bed = combine_list_dfs_to_write(umr_beds, chs)
    # write the results
    df_chs_write_umr_all.to_csv('{}.umr.txt'.format(args.outdir + args.suffix),
                                sep='\t', index=False)
    df_chs_write_umr_bed.to_csv('{}.umr.bed'.format(args.outdir + args.suffix),
                                sep='\t', index=False,
                                header=['# chr', 'chromStart', 'chromEnd', 'name',
                                        'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'])

    # [SUB MODE] save the results
    if args.command == 'UMR-DMR':
      disp('Saving the DMR results')
      df_chs_write_dmr_all = combine_list_dfs_to_write(dmr_all, chs)
      df_chs_write_dmr_bed = combine_list_dfs_to_write(dmr_beds, chs)
      df_chs_write_dmr_all.to_csv('{}.dmr.txt'.format(args.outdir + args.suffix),
                                  sep='\t', index=False)
      df_chs_write_dmr_bed.to_csv('{}.dmr.bed'.format(args.outdir + args.suffix),
                                  sep='\t', index=False,
                                  header=['# chr', 'chromStart', 'chromEnd', 'name',
                                          'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'])
  # [DMR MODE]
  if args.command == 'DMR':
    # raise NotImplementedError('Not implemented yet! Please use UMR-DMR instead!')
    # support multiple umr beds input
    for ibed in args.bed:
      ibed_name = ibed.split('/')[-1].split('.bed')[0]
      names = ['chr', 'pos_start', 'pos_end', 'name', 'score',
      'strand', 'thickStart', 'thickEnd', 'itemRgb']
      df_umr = pd.read_csv(ibed, sep='\t')
      df_umr.columns = names
      chs = df_umr.chr.unique()
      dmr_beds = {}
      for ch in chs:
        # preprocessing the input methylation files
        disp('Processing on {}'.format(ch))
        chr_subs = []
        for i in range(n_files):
          meth_file = meth_files[i]
          chr_sub = meth_file[meth_file.chr == ch]
          chr_sub = merge_strand(chr_sub)
          # downsample the reads proportionally on the sites with high duplicate reads
          chr_sub.loc[:,'methy'], chr_sub.loc[:,'unmethy'] = avoid_overflow_scaling(chr_sub.loc[:,'methy'].values, chr_sub.loc[:,'unmethy'].values)
          chr_subs.append(chr_sub)
        # merge all the data in to a big data frame
        chr_merge = merge_methy_profiles(chr_subs)
        chr_merge = geometric_normalize(chr_merge)

        # preprocess the input umr bed file for each chrom
        df_umr_ch = df_umr[df_umr.chr==ch]
        df_dmr_ch = call_dmr_multi(df_umr_ch, chr_merge, args)
        df_dmr_ch.score = df_dmr_ch.dmr_fdr
        df_dmr_ch.drop(['dmr_fdr'], axis=1, inplace=True)
        dmr_beds[ch] = df_dmr_ch

      # after calling dmr for each chroms; save the dmrs for each ibed file
      # change the column index for writing the bed files
      disp('Saving the DMR results')
      df_chs_write_dmr_bed = combine_list_dfs_to_write(dmr_beds, chs)
      #df_chs_write_dmr_bed.iloc[:,[0,6,7,4,3,5,6,7,8]].to_csv('{}.dmr.bed'.format(args.outdir + ibed_name),sep='\t', index=False,header=['# chr', 'chromStart', 'chromEnd', 'name','score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'])
      df_chs_write_dmr_bed.to_csv('{}.dmr.bed'.format(args.outdir + ibed_name),sep='\t', index=False)

if __name__ == '__main__':
  args = process_args()
  main(args)
