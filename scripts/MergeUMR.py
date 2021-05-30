#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, sys, time, argparse
import pandas as pd
import numpy as np

def disp(text):
    print('[{}] {}'.format(time.asctime(), text), file=sys.stderr)

def Load_UMR(ifile,methratio_max=0.1,FS="\t"):
    disp('Loading UMR: {}'.format(ifile))
    df = pd.read_csv(ifile, sep=FS, header=0)
    df = df.loc[df['mean_methy_ratio']<=methratio_max, ['chr','pos_start','pos_end']]
    df.sort_values(['chr','pos_start','pos_end'], inplace=True, ascending=True)
    #df['prefix'] = pd.Series([prefix] * df.shape[0],index=df.index)
    return df

def read_methy_files(ifile, cols=[0,1,2,6,7]):
    names = ['chr', 'pos', 'strand', 'methy', 'total']
    disp('Loading MethRatio: {}'.format(ifile))
    meth_file = pd.read_csv(ifile, sep='\t', header=0, usecols=cols, names=names, compression='infer')
    meth_file.index = meth_file['pos']
    meth_file.drop(['pos'], axis=1, inplace=True)
    return meth_file

def merge_strand_each_chr(df):
    df_p = df[df['strand']=='+']
    df_n = df[df['strand']=='-']
    df_n.index =  df_n.index.values - 1
    merge_index = np.sort(np.unique(np.append(df_n.index.values, df_p.index.values)))
    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),2]), index=merge_index)
    df_merge.loc[df_p.index,:] = df_merge.loc[df_p.index,:] + df_p.loc[:,['methy','total']].values
    df_merge.loc[df_n.index,:] = df_merge.loc[df_n.index,:] + df_n.loc[:,['methy','total']].values
    df_merge.columns = ['methy','total']
    df_merge = df_merge.loc[0:,:] # remove the minus index pos -1
    return df_merge

def merge_strand(df):
    chs = ['chr' + str(i) for i in range(1, 23)]  + ['chrX','chrY']
    df_merge = pd.DataFrame()
    for ch in chs:
        chr_sub = df[df["chr"] == ch]
        if chr_sub.shape[0] > 0:
            #disp('Merging MethRatio on {}'.format(ch))
            chr_sub = merge_strand_each_chr(chr_sub)
            chr_sub['chr']=pd.Series([ch] * chr_sub.shape[0], index=chr_sub.index)
            df_merge=df_merge.append(chr_sub)
    return df_merge

def Region_Meth_Ratio(methratio_sub,start=0,end=0):
    #methratio_sub: methratio df of one chrom
    region_meth=methratio_sub.loc[start:(end+1),:]
    count_C=region_meth.shape[0]
    methy_C=region_meth["methy"].sum()
    total_C=region_meth["total"].sum()
    if total_C==0:
        region_methratio=np.nan
    else:
        region_methratio=methy_C*1.0/total_C
    return region_methratio,count_C

def Merge_Adjacent_UMR(umr,methratio,adjacent_max=500,methratio_max=0.1,prefix="out"):
    umr_merge=pd.DataFrame()
    chs = ['chr' + str(i) for i in range(1, 23)]  + ['chrX','chrY']
    for ch in chs:
        umr_a=umr[umr['chr']==ch]
        if umr_a.shape[0] > 0:
            #disp('Merging Adjacent UMR on {}'.format(ch))
            methratio_sub=methratio[methratio['chr']==ch]
            start0=0;end0=0;
            for row in umr_a.iterrows():
                start1=row[1]['pos_start'];
                end1=row[1]['pos_end'];
                if end0==0:
                    start0=start1;end0=end1;
                else:
                    if start1-end0 <= adjacent_max:
                        mratio1,count1=Region_Meth_Ratio(methratio_sub=methratio_sub,start=start0,end=end1)
                        if mratio1 <= methratio_max:
                            end0=end1;
                        else:
                            umr_merge=umr_merge.append(pd.DataFrame([[ch,int(start0),int(end0)]]))
                            start0=start1;end0=end1;
                    else:
                        umr_merge=umr_merge.append(pd.DataFrame([[ch,int(start0),int(end0)]]))
                        start0=start1;end0=end1;
            umr_merge=umr_merge.append(pd.DataFrame([[ch,int(start0),int(end0)]]))
    umr_merge.columns = ['chr','pos_start','pos_end']
    disp('Output Merged UMR: n={}'.format(umr_merge.shape[0]))
    umr_merge.to_csv(prefix+"_Merged_UMR.bed",sep="\t",columns=['chr','pos_start','pos_end'],header=False,index=False)
    #return umr_merge

ScriptName=re.split('/', sys.argv[0])[-1]
cmds={"MergeUMR":"Merge Adjacent UMRs if they meet specified threshold."}

def printHelp():
    print("For help information of this script, try:\n")
    print("  python "+ScriptName+" <Function> -h\n")
    print("Availible Functions:\n")
    for i in cmds.keys():
        print('  '+i+'\t'+cmds[i]+'\n')

def main(cmd=''):
    # check whether the requested function is included
    if cmd not in cmds.keys():
        print("Function ",cmd," not found, please see the help below:\n")
        printHelp()
        #sys.exit()
        return False

    # parse parameters
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage='python '+ScriptName+' '+cmd+' <positional arguments> [optional arguments]\n',\
                                     description=cmds[cmd]+"\n",\
                                     epilog='')
    # designate the first positional argument as ' '
    parser.add_argument(' ',default=None,help='')

    # parse input parameters of each function
    if cmd=='MergeUMR':
        parser.add_argument('UMR',default=None,\
                            help="UMRs called by PopCanyon")
        parser.add_argument('MethRatio',default=None,\
                            help="MethRatio file from BSMAP which was used to call the UMR")
        parser.add_argument('-r', '--methratio_max',dest="methratio_max",metavar='',default=0.1,type=float,\
                            help="Max methratio for UMR detection and Adjacent UMR merging")
        parser.add_argument('-a', '--adjacent_max',dest="adjacent_max",metavar='',default=500,type=int,\
                            help="Max distance of merged adjacent UMR")
        parser.add_argument('-p', '--prefix',dest="OUT_prefix",metavar='',default='out',\
                            help="Prefix of output file")

    # print help information for the requested function
    if '-h' in sys.argv or '--help' in sys.argv:
        parser.print_help()
        return False
    else:
        if len(sys.argv)==2:
            print("Too few arguments. Try:\npython "+ScriptName+" "+cmd+" -h")
            return False

    # save all paramter values in args
    args=parser.parse_args()

    # run function
    if cmd=='MergeUMR':
        disp("MergeUMR Started")
        UMR_df=Load_UMR(ifile=args.UMR,methratio_max=args.methratio_max)
        MethRatio_df=read_methy_files(ifile=args.MethRatio, cols=[0,1,2,6,7])
        MethRatio_df=merge_strand(df=MethRatio_df)
        UMR_df=Merge_Adjacent_UMR(umr=UMR_df,methratio=MethRatio_df,\
                                  adjacent_max=args.adjacent_max,\
                                  methratio_max=args.methratio_max,\
                                  prefix=args.OUT_prefix)
        disp("MergeUMR Finished")

if len(sys.argv)>1:
    main(cmd=sys.argv[1])
else:
    printHelp()
