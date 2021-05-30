import pandas as pd
import numpy as np
import os
import sys
import time
import scipy.stats as ss
from scipy.special import digamma, gammaln
from scipy.optimize import minimize
try:
    from rpy2 import robjects
    robject_failure = False
except:
    robject_failure = True

################################### Below Functions for HMM Beta-Binomial ###################################

#read all the methy files for each chroms to reduce the memory consumption
# return pandas dataframe for each chrom
def read_methy_files_old(files, ch, cols=[0,1,2,6,7]):
    names = ['chr', 'pos', 'strand', 'methy', 'total']
    meth_files = []
    for ifile in files:
        with open(ifile, 'r') as f:
            meth_file_dict = {name:[] for name in names}
            header = f.readline().rstrip().split('\t')
            # if the first line is header then omit it
            try:
                pos = int(header[1])
                if header[0] == ch:
                    meth_file_dict['chr'].append(header[cols[0]])
                    meth_file_dict['pos'].append(int(header[cols[1]]))
                    meth_file_dict['strand'].append(header[cols[2]])
                    meth_file_dict['methy'].append(int(header[cols[3]]))
                    meth_file_dict['total'].append(int(header[cols[4]]))
            except: pass

            for line in f:
                eles = line.rstrip().split('\t')
                if eles[0] == ch:
                    meth_file_dict['chr'].append(eles[cols[0]])
                    meth_file_dict['pos'].append(int(eles[cols[1]]))
                    meth_file_dict['strand'].append(eles[cols[2]])
                    meth_file_dict['methy'].append(int(eles[cols[3]]))
                    meth_file_dict['total'].append(int(eles[cols[4]]))
            meth_file = pd.DataFrame(meth_file_dict)
            meth_file = meth_file.reindex(columns = names)
            meth_file.index = meth_file['pos']
            meth_file.drop(['pos'], axis=1, inplace=True)
        meth_files.append(meth_file)
    return meth_files
############### the above function save memory but slow; will be removed in future
#read all the methy files for all chroms
# index the files with the pos
def read_methy_files(files, cols=[0,1,2,6,7]):
    names = ['chr', 'pos', 'strand', 'methy', 'total']
    meth_files = []
    for ifile in files:
        disp('Reading file: {}'.format(ifile))
        meth_file = pd.read_csv(ifile, sep='\t', header=1, usecols=cols, names=names, compression='infer')
        meth_file.index = meth_file['pos']
        meth_file.drop(['pos'], axis=1, inplace=True)
        meth_files.append(meth_file)

    return meth_files


# return Ture if it is nan
def is_nan(x):
    return (x is np.nan or x != x)


# avoid the overflow when coming across large number of reads
def avoid_overflow_scaling(x_data, y_data, overflow_val=200):
    overflow_id = np.argwhere(x_data > overflow_val)
    while len(overflow_id) > 0:
        # print(overflow_id, 'x')
        x_data[overflow_id[:,0]] = x_data[overflow_id[:,0]] / 2
        y_data[overflow_id[:,0]] = y_data[overflow_id[:,0]] / 2
        overflow_id = np.argwhere(x_data > overflow_val)

    overflow_id = np.argwhere(y_data > overflow_val)
    while len(overflow_id) > 0:
        # print(overflow_id, 'y')
        y_data[overflow_id[:,0]] = y_data[overflow_id[:,0]] / 2
        x_data[overflow_id[:,0]] = x_data[overflow_id[:,0]] / 2
        overflow_id = np.argwhere(y_data > overflow_val)

    x_data = np.ceil(x_data)
    y_data = np.ceil(y_data)

    return x_data, y_data


## add rep to the likelihood
## betabinomial (bb) liklihood for kth component
def compute_minus_l_rep_mix(x,post,x_data,y_data):
    M = x_data.shape[1]
    l = M*post.sum()*( gammaln(sum(x)) - np.sum(gammaln(x)) )\
     + post.dot(gammaln(x_data+x[0]).sum(axis=1))\
     + post.dot(gammaln(y_data+x[1]).sum(axis=1))\
     - post.dot(gammaln(x_data+y_data+np.sum(x)).sum(axis=1))
    return -l


## bb derivative of likelihood for kth component
def compute_J(x,post,x_data,y_data):
    M = x_data.shape[1]
    J = M*post.sum()*digamma(sum(x)) \
    - M*post.sum()*digamma(x) \
    - post.dot(digamma(x_data+y_data+np.sum(x)).sum(axis=1)) \
    + post.dot(np.concatenate([digamma(x_data+x[0]).sum(axis=1).reshape(-1,1),
        digamma(y_data+x[1]).sum(axis=1).reshape(-1,1)],axis=1))
    return -J


# the basic version for solving the hmm betabinomial with replicates
def hmm_bb_mix3_rep_solver(x_data,y_data,Nit=10):
    para = np.array([[5,1],[1,5],[5,5]],dtype=np.float) ##initializer
    alpha = para[:, 0]
    beta = para[:, 1]
    K, K_para = para.shape
    x_data += 1
    y_data += 1
    N, M = x_data.shape
    trans = np.array([[1.0/K]*K for i in range(K)])
    forward = np.zeros((N,K),dtype=np.float64)
    backward = np.zeros((N,K),dtype=np.float64)
    dens = np.zeros((N,K),dtype=np.float64)
    scale = np.array([1] + [0] * (N - 1),dtype=np.float64)
    bnds = tuple((0,1000) for _ in range(K_para))
    postprob = np.zeros((N,K),dtype=np.float)

    # x_data, y_data = avoid_overflow_scaling(x_data, y_data)
    ## M-step for optimizing the dirichlet parameters
    for iteration in range(Nit):
        ## density function with replicates
        dens = np.exp(M*gammaln(alpha+beta) - M*gammaln(alpha) - M*gammaln(beta) +
        gammaln(x_data.reshape(-1,M,1)+alpha).sum(axis=1) + gammaln(y_data.reshape(-1,M,1)+beta).sum(axis=1) - gammaln((x_data+y_data).reshape(-1,M,1)+alpha+beta).sum(axis=1))

        ## to avoid the overflow for some extreme large data
        ## in most cases will not come across such situation
        overflow_id = np.argwhere(dens.sum(1) == 0).ravel()
        while len(overflow_id) > 1:
            print('overflow')
            x_data[overflow_id] = np.ceil(x_data[overflow_id] / 10)
            y_data[overflow_id] = np.ceil(y_data[overflow_id] / 10)
            dens = np.exp(M*gammaln(alpha+beta) - M*gammaln(alpha) - M*gammaln(beta) +
            gammaln(x_data.reshape(-1,M,1)+alpha).sum(axis=1) + gammaln(y_data.reshape(-1,M,1)+beta).sum(axis=1) - gammaln((x_data+y_data).reshape(-1,M,1)+alpha+beta).sum(axis=1))
            overflow_id = np.argwhere(dens.sum(1) == 0).ravel()

        ## E-step for postprob
        ## forward recursion
        forward[0,:] = dens[0,:] / K
        for t in range(1,N):
            forward[t,:] = forward[t-1,:].dot(trans) * dens[t,:]
            scale[t] = forward[t,:].sum()
            forward[t,:] = forward[t,:] / scale[t]

        print(np.log(sum(forward[N-1,:])) + sum(np.log(scale)))

        ## backward recursion
        backward[N-1,:] = np.array([1]*K)
        for t in reversed(range(N-1)):
            backward[t,:] = (backward[t+1,:] * dens[t+1,:]).dot(trans.T)
            backward[t,:] /= scale[t]

        ## transition reparametrization
        trans = trans * (forward[:N-1,:].T).dot(dens[1:,:] * backward[1:,:])
        trans = trans / (trans.sum(axis=1).reshape(-1,1).dot(np.ones([1,K])))

        ## postprob reparametrization
        postprob = forward * backward
        postprob = postprob / (postprob.sum(axis=1).reshape(-1,1).dot(np.ones([1,K])))
        wghts = postprob.sum(axis=0)
        wghts /= wghts.sum()

        # compute the parameters
        for j in range(K):
            x0 = para[j]
            res = minimize(compute_minus_l_rep_mix, x0, args=(postprob[:,j],x_data,y_data,), method='SLSQP', jac=compute_J,hess=None, bounds=bnds, tol=None, callback=None, options={'disp': False, 'eps': 1.4901161193847656e-08, 'maxiter': 5 ,'ftol':0.01})
            if not any(is_nan(res.x)):
                para[j] = res.x
            else:
                break

    return (para,postprob,wghts)
# TODO: Above function will be removed in future
################################### Above Functions for HMM Beta-Binomial ##########################

################################### Below Functions for Diff analysis ###############################
# bb solver for single component
def bb_rep_solver(x_data, y_data, Nit=15):
    x0 = np.array([5,5])
    post = np.ones(shape=(x_data.shape[0],))
    bnds = tuple((0,1000) for _ in range(len(x0)))
    res_fit = minimize(compute_minus_l_rep_mix, x0, args=(post,x_data,y_data,), method='SLSQP', jac=compute_J,hess=None, bounds=bnds, tol=None, callback=None, options={'disp': False, 'eps': 1.4901161193847656e-08, 'maxiter': Nit ,'ftol':0.01})
    return res_fit

# bb likelihood ratio test;
# umr_counts have the format: x1 y1 x2 y2 x3 y3
# treat_num denote how many treats samples with treat samples first
def bb_likelihood_ratio_test(umr_counts, labels):
    treat_id = [i for i, label in enumerate(labels) if label=='treat']
    ctrl_id = [i for i, label in enumerate(labels) if label=='ctrl']
    if len(treat_id) < 1:
        raise ValueError('Could not find treat samples!')
    if len(ctrl_id) < 1:
        raise ValueError('Could not find ctrl samples!')
    counts = umr_counts.values.copy()
    counts += .2 # to avoid dividing by zero
    x_data = counts[:,range(0,counts.shape[1],2)]
    y_data = counts[:,range(1,counts.shape[1],2)]
    res_treat_ctrl = bb_rep_solver(x_data, y_data,30)
    res_treat = bb_rep_solver(x_data[:,treat_id], y_data[:,treat_id])
    res_ctrl = bb_rep_solver(x_data[:,ctrl_id], y_data[:,ctrl_id])
    tst = (res_treat_ctrl.fun - res_treat.fun - res_ctrl.fun) * 2
    # p_val = 1 - robjects.r['pchisq'](tst, 2)[0] # not significant faster than scipy
    p_val = 1 - ss.chi2.cdf(tst,df=2)
    # calcualte the differential methylation values
    treat_methy_ratio = (x_data[:,treat_id] / (x_data[:,treat_id] + y_data[:,treat_id]))
    ctrl_methy_ratio = (x_data[:,ctrl_id] / (x_data[:,ctrl_id] + y_data[:,ctrl_id]))
    diff_methy_ratio = treat_methy_ratio - ctrl_methy_ratio
    diff_mean_methy_ratio = np.mean(diff_methy_ratio)
    diff_var_methy_ratio = np.var(diff_methy_ratio)

    return p_val, diff_mean_methy_ratio, diff_var_methy_ratio


#################################### Above Functions for Diff analysis ##############################

################################### Below Functions for Utilities ###################################
# display the intermediate steps and results
def disp(txt):
    print('[@{}] {}'.format(time.asctime(), txt), file=sys.stderr)

# cacluate the fdr using pvals
def pval2fdr(pval):
    pval = robjects.FloatVector(pval)
    fdr = robjects.r['p.adjust'](pval,'fdr')
    return np.array(fdr)


## merge positive strand and negative strand reads together
#def merge_strand(df):
#    df_p = df[df['strand']=='+']
#    df_n = df[df['strand']=='-']
#    df_n.index =  df_n.index.values - 1
#    merge_index = np.sort(np.unique(np.append(df_n.index.values, df_p.index.values)))
#    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),2]), index=merge_index)
#    df_merge.loc[df_p.index,:] = df_merge.loc[df_p.index,:] + df_p.loc[:,['methy','total']].values
#    df_merge.loc[df_n.index,:] = df_merge.loc[df_n.index,:] + df_n.loc[:,['methy','total']].values
#    df_merge.columns = ['methy','total']
#    df_merge['unmethy'] = df_merge['total'] - df_merge['methy']
#    df_merge = df_merge.reindex(columns=['methy','unmethy','total'])
#    return df_merge

# TODO: if no bugs in new merge_strand, then the above one will be removed
# merge positive strand and negative strand reads together
def merge_strand(df):
    df_p = df[df['strand']=='+']
    df_n = df[df['strand']=='-']
    df_n.index =  df_n.index.values - 1
    merge_index = np.sort(np.unique(np.append(df_n.index.values, df_p.index.values)))
    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),2]), index=merge_index)
    df_merge.loc[df_p.index,:] = df_merge.loc[df_p.index,:] + df_p.loc[:,['methy','total']].values
    df_merge.loc[df_n.index,:] = df_merge.loc[df_n.index,:] + df_n.loc[:,['methy','total']].values
    df_merge.columns = ['methy','total']
    df_merge['unmethy'] = df_merge['total'] - df_merge['methy']
    df_merge = df_merge.reindex(columns=['methy','unmethy'])
    df_merge = df_merge.loc[0:,:] # remove the minus index pos -1
    return df_merge


# merge list of methylation data together
# L: list of the methylation profile with columns named as
# methy, total; each methylation profile should be indexed according to the pos; .loc include the end position for slicing
def merge_methy_profiles(L):
    all_index = np.array([],dtype=np.int)
    for l in L:
        all_index = np.append(all_index,l.index.values)
    # K: denote the number of cols in each methy file; Here only includes methy and unmethy
    K = 2
    merge_index =  np.sort(np.unique(all_index))
    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),K*len(L)]), index=merge_index)
    for i, l in enumerate(L):
        df_merge.loc[l.index,i*K:(i*K)+(K-1)] = l.loc[:,['methy','unmethy']].values
    return df_merge


# column order for df should be always methy, unmethy, methy, unmethy
def geometric_normalize(df):
    total_counts = []
    for i in range(0,df.shape[1],2):
        total_counts.append(np.sum(df.iloc[:,i:i+2].sum(1)))
    total_counts = np.array(total_counts)
    gmean = np.power(np.prod(total_counts), 1.0/len(total_counts))
    normalizer = gmean / total_counts
    for i in range(len(normalizer)):
        df.iloc[:,i*2:i*2+2] =  normalizer[i] * df.iloc[:,i*2:i*2+2]
    return np.round(df)


def extract_consecutive_region(vec, max_interval=2):
    vec_interval = vec[1:] - vec[:-1]
    vec_mark = np.where(vec_interval > max_interval)[0] + 1
    return np.insert(vec_mark, [0], [0])


def pbinom_slow(vec, r):
    return ss.binom_test(vec[0], vec[1], r, 'less')


# Use vectorized version of binomial test
def pbinom(vec, vec_total, r):
    vec = robjects.IntVector(vec)
    vec_total = robjects.IntVector(vec_total)
    p_val = robjects.r['pbinom'](vec, vec_total, r)
    return np.array(p_val)


# check whether the directory exits if not create one
def check_and_create_dir(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)


# according to the viterbi returned class and fdr, get one line of the region from pos start to pos end
# use the peak region to re-do the p-val test for the whole region
def get_region_from_class_label(df_comb, args):
    df_comb_res = df_comb[df_comb['umr_class'].values > 0]
    id2index = {i:index for i,index in enumerate(df_comb_res.index)}
    id2pos = {i:pos for i,pos in enumerate(df_comb_res.pos)}
    # stores the methy ratio for each sample
    all_methy_ratio = df_comb_res.methy.values / (df_comb_res.methy.values + df_comb_res.unmethy.values + 0.1)
    id2methy_ratio = {i:all_methy_ratio[i] for i in range(all_methy_ratio.shape[0])}
    id2methy_counts = {i:counts for i, counts in enumerate(df_comb_res.methy_all)}
    id2unmethy_counts = {i:counts for i, counts in enumerate(df_comb_res.unmethy_all)}
    region_index_bb = extract_consecutive_region(df_comb_res.index,max_interval=2)
    umr_dict = {'pos_start':[],'pos_end':[],\
                'index_start':[], 'index_end':[], 'num_methy_sites':[], \
                'region_width':[], 'mean_methy_ratio':[], 'variance_methy_ratio':[],\
                'total_methy_counts':[], 'total_unmethy_counts':[]}
    for i in range(len(region_index_bb)-1):
        if region_index_bb[i + 1] -  region_index_bb[i] >= args.min_sites:
            # region_index_bb[i+1]: next region start loci
            umr_dict['pos_start'].append(id2pos[region_index_bb[i]])
            umr_dict['pos_end'].append(id2pos[region_index_bb[i+1]-1])
            umr_dict['index_start'].append(id2index[region_index_bb[i]])
            umr_dict['index_end'].append(id2index[region_index_bb[i+1]-1])
            umr_dict['num_methy_sites'].append(id2index[region_index_bb[i+1]-1] - id2index[region_index_bb[i]] + 1)
            umr_dict['region_width'].append(id2pos[region_index_bb[i+1]-1] - id2pos[region_index_bb[i]] + 1)
            # add more detailed stats on the anaylsis
            methy_counts = []
            unmethy_counts = []
            methy_ratio = []
            for j in range(region_index_bb[i], region_index_bb[i+1]):
                methy_counts.append(id2methy_counts[j])
                unmethy_counts.append(id2unmethy_counts[j])
                methy_ratio.append(id2methy_ratio[j])

            umr_dict['total_methy_counts'].append(np.array(methy_counts,dtype=int).sum())
            umr_dict['total_unmethy_counts'].append(np.array(unmethy_counts,dtype=int).sum())
            umr_dict['mean_methy_ratio'].append(np.array(methy_ratio).ravel().mean())
            umr_dict['variance_methy_ratio'].append(np.array(methy_ratio).ravel().var())


    return pd.DataFrame(umr_dict)


# according to the index and fdr, get one line of the region from pos start to pos end
def get_region_from_index(df_comb_res, which_fdr, args):
    # extract the region for consecutive significant FDR
    id2index = {i:index for i,index in enumerate(df_comb_res.index)}
    id2pos = {i:pos for i,pos in enumerate(df_comb_res.pos)}
    id2fdr = {i:fdr for i,fdr in enumerate(df_comb_res.loc[:,which_fdr])}
    # stores the methy ratio for each sample
    all_methy_ratio = df_comb_res.methy.values / (df_comb_res.methy.values + df_comb_res.unmethy.values)
    id2methy_ratio = {i:all_methy_ratio[i] for i in range(all_methy_ratio.shape[0])}
    id2methy_counts = {i:counts for i, counts in enumerate(df_comb_res.methy_all)}
    id2unmethy_counts = {i:counts for i, counts in enumerate(df_comb_res.unmethy_all)}
    region_index_bb = extract_consecutive_region(df_comb_res.index.values)
    umr_dict = {'pos_start':[],'pos_end':[],'umr_fdr':[],\
                'index_start':[], 'index_end':[], 'num_methy_sites':[], \
                'region_width':[], 'mean_methy_ratio':[], 'variance_methy_ratio':[],\
                'total_methy_counts':[], 'total_unmethy_counts':[]}
    for i in range(len(region_index_bb)-1):
        if region_index_bb[i + 1] -  region_index_bb[i] >= args.min_sites:
            # region_index_bb[i+1]: next region start loci
            umr_dict['pos_start'].append(id2pos[region_index_bb[i]])
            umr_dict['pos_end'].append(id2pos[region_index_bb[i+1]-1])
            umr_dict['index_start'].append(id2index[region_index_bb[i]])
            umr_dict['index_end'].append(id2index[region_index_bb[i+1]-1])
            umr_dict['num_methy_sites'].append(id2index[region_index_bb[i+1]-1] - id2index[region_index_bb[i]] + 1)
            umr_dict['region_width'].append(id2pos[region_index_bb[i+1]-1] - id2pos[region_index_bb[i]] + 1)
            # add more detailed stats on the anaylsis
            methy_counts = []
            unmethy_counts = []
            fdrs = []
            methy_ratio = []
            for j in range(region_index_bb[i], region_index_bb[i+1]):
                methy_counts.append(id2methy_counts[j])
                unmethy_counts.append(id2unmethy_counts[j])
                fdrs.append(id2fdr[j])
                methy_ratio.append(id2methy_ratio[j])
            umr_dict['total_methy_counts'].append(np.array(methy_counts,dtype=int).sum())
            umr_dict['total_unmethy_counts'].append(np.array(unmethy_counts,dtype=int).sum())
            umr_dict['umr_fdr'].append(np.array(fdrs).mean())
            umr_dict['mean_methy_ratio'].append(np.array(methy_ratio).ravel().mean())
            umr_dict['variance_methy_ratio'].append(np.array(methy_ratio).ravel().var())


    return pd.DataFrame(umr_dict)


########################## slow version and the original version
def get_region_from_index_slow(df_comb_res, which_fdr, args):
    # extract the region for consecutive significant FDR
    region_index_bb = extract_consecutive_region(df_comb_res.index.values)
    umr_dict = {'pos_start':[],'pos_end':[],'umr_fdr':[],\
                'index_start':[], 'index_end':[]}
    for i in range(len(region_index_bb)-1):
        if region_index_bb[i + 1] -  region_index_bb[i] >= args.min_cpg_sites:
            # region_index_bb[i+1]: next region start loci
            umr_dict['pos_start'].append(df_comb_res.iloc[region_index_bb[i]].pos.astype('int'))
            umr_dict['pos_end'].append(df_comb_res.iloc[region_index_bb[i+1]-1].pos.astype('int'))
            umr_dict['fdr'].append(df_comb_res.iloc[region_index_bb[i]:region_index_bb[i+1]-1].loc[:,which_fdr].mean())
            umr_dict['index_start'].append(df_comb_res.iloc[region_index_bb[i]].name)
            umr_dict['index_end'].append(df_comb_res.iloc[region_index_bb[i+1]-1].name)
    return pd.DataFrame(umr_dict)
############################### will be removed in future


################################ The HMM beta model for the single sample

# call dmrs based on the previously detected umrs
# add the differential mean ratio value and variance in bb_likelihood_ratio function
def call_dmr_multi(df_umr, chr_merge, args):
    umr = df_umr.copy()
    diff_pval = []
    diff_mean_methy_ratio = []
    diff_var_methy_ratio = []
    for start, end in zip(umr.pos_start, umr.pos_end):
        umr_counts = chr_merge.loc[start:end]
        p_val, diff_mean_methy, diff_var_methy = bb_likelihood_ratio_test(umr_counts, args.labels)
        diff_pval.append(p_val)
        diff_mean_methy_ratio.append(diff_mean_methy)
        diff_var_methy_ratio.append(diff_var_methy)

    # store the fdr for the dmr
    umr['dmr_fdr'] = pval2fdr(np.array(diff_pval))
    umr['diff_mean_methy_ratio'] = np.array(diff_mean_methy_ratio)
    umr['diff_var_methy_ratio'] = np.array(diff_var_methy_ratio)
    umr = umr[umr.dmr_fdr <= args.dmr_fdr]

    # write the color
    hyper = umr.diff_mean_methy_ratio > 0.05
    hypo = umr.diff_mean_methy_ratio < -0.05
    umr.loc[hyper,'itemRgb'] = '244,107,66'
    umr.loc[hypo, 'itemRgb'] = '66,244,179'

    return umr

# given the df_umr, which includes methylation start, end and raw data, compute
# the statistics for each detected umr
# return: averaged methylation ratio and its variances;
def get_region_stats(df_umr, df_comb_res_umr):
    umr_dict = {'averaged_methy_ratio':[],'variance_methy_ratio':[],'total_methy_reads':[],\
    'total_unmethy_reads':[], 'number_cpg_sites':[]}
    for i in range(df_umr.shape[0]):
        df_tmp = df_comb_res_umr.loc[df_umr.index_start[i]:df_umr.index_end[i]]
        methy_ratio = df_tmp.methy.values / (df_tmp.methy.values + df_tmp.unmethy.values)
        umr_dict['averaged_methy_ratio'].append(np.mean(methy_ratio.ravel()))
        umr_dict['variance_methy_ratio'].append(np.var(methy_ratio.ravel()))
        umr_dict['total_methy_reads'].append(np.sum(df_tmp.methy_all.values))
        umr_dict['total_unmethy_reads'].append(np.sum(df_tmp.unmethy_all.values))
        umr_dict['number_cpg_sites'].append(df_umr.index_end[i]-df_umr.index_start[i]+1)
    umr_stats = pd.DataFrame(umr_dict)
    umr_stats.reindex(columns=['number_cpg_sites', 'averaged_methy_ratio', 'variance_methy_ratio', 'total_methy_reads', 'total_unmethy_reads'])
    return umr_ststs

####### above get_region_stats is integrated into the get_region_from_index and is much faster
####### will be removed in future

# use viterbi to decode the states
def viterbi_decoder(para, x_data, y_data, trans):
    alpha = para[:, 0]
    beta = para[:, 1]
    K = para.shape[1]
    N, M = x_data.shape
    dens = np.exp(M*gammaln(alpha+beta) - M*gammaln(alpha) - M*gammaln(beta) +
        gammaln(x_data.reshape(-1,M,1)+alpha).sum(axis=1) + gammaln(y_data.reshape(-1,M,1)+beta).sum(axis=1) - gammaln((x_data+y_data).reshape(-1,M,1)+alpha+beta).sum(axis=1))
    logdens = np.log(dens + 2.2e-300) # 2.2e-300: very small double number to avoid overflow
    trans = np.log(trans + 2.2e-300)
    # dynamic programming
    plogl = np.zeros(shape=(N,K))
    backtr = np.zeros(shape=(N-1,K))
    plogl[0,:] = logdens[0,:] - np.log(K)

    for t in range(1, N):
        tmp = plogl[t-1].reshape(-1,1).dot(np.ones(shape=(1,K))) + trans
        plogl[t,:] = tmp.max(axis=0)
        backtr[t-1,:] = np.argmax(tmp, axis=0)
        plogl[t,:] = plogl[t,:] + logdens[t,:]

    clss = np.zeros(shape=(N,), dtype=int)
    clss[N-1] = np.argmax(plogl[N-1])
    # find the hidden status
    for t in range(N-2,-1,-1):
        clss[t] = backtr[t,clss[t+1]]

    return clss


# solve the bb mixture with HMM for K components
# supporting to solve 2 or 3 mixture components by default
def hmm_bb_mix_rep_solver_K(x_data,y_data,Nit=10,K=2):
    if K == 2:
        para = np.array([[5,1],[1,5]],dtype=np.float) ##initializer
    if K == 3:
        para = np.array([[5,1],[1,5],[15,15]],dtype=np.float)
    alpha = para[:, 0]
    beta = para[:, 1]
    K_para = para.shape[1]
    x_data = x_data.copy()
    y_data = y_data.copy()
    x_data += 1
    y_data += 1
    N, M = x_data.shape
    trans = np.array([[1.0/K]*K for i in range(K)])
    forward = np.zeros((N,K),dtype=np.float64)
    backward = np.zeros((N,K),dtype=np.float64)
    dens = np.zeros((N,K),dtype=np.float64)
    scale = np.array([1] + [0] * (N - 1),dtype=np.float64)
    bnds = tuple((0,1000) for _ in range(K_para))
    postprob = np.zeros((N,K),dtype=np.float)

    ## M-step for optimizing the dirichlet parameters
    for iteration in range(Nit):
        ## density function with replicates
        dens = np.exp(M*gammaln(alpha+beta) - M*gammaln(alpha) - M*gammaln(beta) +
        gammaln(x_data.reshape(-1,M,1)+alpha).sum(axis=1) + gammaln(y_data.reshape(-1,M,1)+beta).sum(axis=1) - gammaln((x_data+y_data).reshape(-1,M,1)+alpha+beta).sum(axis=1))

        ## forward recursion
        forward[0,:] = dens[0,:] / K
        for t in range(1,N):
            forward[t,:] = forward[t-1,:].dot(trans) * dens[t,:]
            scale[t] = forward[t,:].sum()
            forward[t,:] = forward[t,:] / scale[t]

        # print(np.log(sum(forward[N-1,:])) + sum(np.log(scale)))

        backward = np.zeros((N,K),dtype=np.float64)
        ## backward recursion
        backward[N-1,:] = np.array([1]*K)
        for t in reversed(range(N-1)):
            backward[t,:] = (backward[t+1,:] * dens[t+1,:]).dot(trans.T)
            backward[t,:] /= scale[t]
            if not all(np.isfinite(backward[t,:])):
                print(t); break;

        ## transition reparametrization
        trans = trans * (forward[:N-1,:].T).dot(dens[1:,:] * backward[1:,:])
        trans = trans / (trans.sum(axis=1).reshape(-1,1).dot(np.ones([1,K])))

        ## postprob reparametrization
        postprob = forward * backward
        postprob = postprob / (postprob.sum(axis=1).reshape(-1,1).dot(np.ones([1,K])))
        wghts = postprob.sum(axis=0)
        wghts /= wghts.sum()

        # compute the parameters
        for j in range(K):
            x0 = para[j]
            res = minimize(compute_minus_l_rep_mix, x0, args=(postprob[:,j],x_data,y_data,), method='SLSQP', jac=compute_J,hess=None, bounds=bnds, tol=None, callback=None, options={'disp': False, 'eps': 1.4901161193847656e-08, 'maxiter': 15 ,'ftol':0.01})
            if not any(is_nan(res.x)):
                para[j] = res.x
            else:
                break

    return (para,postprob,wghts,trans)


# HMM and Beta-binomial computation framework
# df_sub: all the methy unmethy info with pvals meeting the threshold
# region_index: consecutive region indices
# df_backup: original methy unmethy counts without sliding
# This version could handle input samples with / without replicates
def BB_HMM_compute_multi(df_sub, df_backup, args):
    # modified for wrapping the select the region into a fucntion get_region_from_index
    # extract the hotspot region which is the smallest pval region plus the neighboring region
    id2index = {i:index for i,index in enumerate(df_sub.index)}
    region_index = extract_consecutive_region(df_sub.index.values)
    df_tmp_list = []
    for i in range(len(region_index)-1):
        if region_index[i + 1] -  region_index[i] >= args.min_sites:
            # expand the hotspot region to two directions
            hotspot_index_start = id2index[region_index[i]]
            hotspot_index_end = id2index[region_index[i+1]-1]
            hotspot_len = hotspot_index_end - hotspot_index_start + 1
            index_start = hotspot_index_start - 2*hotspot_len
            index_end = hotspot_index_end + 2*hotspot_len
            df_tmp = df_backup.iloc[index_start:index_end]
            df_tmp_list.append(df_tmp)

    # no umr found in this region
    if len(df_tmp_list) == 0:
        # return empty dataframe
        df_umr = pd.DataFrame([], columns=['chr', 'pos_start', 'pos_end', 'index_start', 'index_end', 'umr_fdr', 'num_methy_sites', 'region_width', 'mean_methy_ratio', 'variance_methy_ratio', 'total_methy_counts', 'total_unmethy_counts'])
        return df_umr

    if len(df_tmp_list) > 1:
        df_comb = df_tmp_list[0].append(df_tmp_list[1:])
    else:
        df_comb = df_tmp_list[0]

    uniq_id = np.unique(df_comb.index, return_index=True)[1]
    # df_comb: stores all the methy unmethy, index, pos info
    df_comb = df_comb.iloc[uniq_id,:]
    print('Computing the umr using Beta-binomial HMM model')
    N = df_comb.shape[0]
    # more than one samples
    x_data = df_comb.methy.values
    y_data = df_comb.unmethy.values
    # there is only one sample in x_data
    # row vector x_data should be converted to column vector
    if len(x_data.shape) == 1:
        x_data = df_comb.methy.values.reshape(N,1)
        y_data = df_comb.unmethy.values.reshape(N,1)

    ab, post, wght, trans = hmm_bb_mix_rep_solver_K(x_data+1, y_data+1, Nit=args.Nit, K=2)
    clss = viterbi_decoder(ab, x_data, y_data, trans)

#    ab, post, wght = hmm_bb_mix3_rep_solver_inter_fix(df_comb.methy.values,df_comb.unmethy.values)
    # according to which component (UMR or IMR), assign the postprob
    # ab: the fitted alpha beta values for the 3 mixture beta-binomial distribution
    methy_stats = ab[:,0] / ab.sum(1)
    umr_ind, imr_ind = np.argsort(methy_stats)[:2]
    df_comb['umr_class'] = np.array(clss == umr_ind).astype(int)
    df_comb['umr_fdr'] = pval2fdr(1 - post[:,umr_ind])
    df_comb['umr_post'] = post[:,umr_ind]
    # 2018-2-16: add the option to compute the rectfied fdr
    if not args.viterbi_disable:
        df_umr = get_region_from_class_label(df_comb, args)
        # df_umr['umr_fdr'] = pval2fdr(calculate_region_pval(df_umr, args))
        df_umr['umr_fdr'] = pval2fdr(calculate_region_pval(df_umr, args))
        df_umr = df_umr[df_umr.umr_fdr < args.umr_fdr]
    else:
        df_comb_res_umr = df_comb[df_comb.umr_fdr <= args.umr_fdr]
        df_umr = get_region_from_index(df_comb_res_umr,'umr_fdr', args)
    return df_umr


# fdr is too small for testing the region, so we normalize the region length
# def calculate_region_pval(df_umr, args):
#     avg_counts = (df_umr.loc[:,'total_methy_counts'].sum() + df_umr.loc[:,'total_unmethy_counts'].sum()) /\
#                     df_umr.loc[:,'region_width'] .sum()
#     rect_methy = df_umr.loc[:,'total_methy_counts'] / df_umr.loc[:,'region_width'] * avg_counts
#     rect_unmethy = df_umr.loc[:,'total_unmethy_counts'] / df_umr.loc[:,'region_width'] * avg_counts
#     pval = pbinom(rect_methy, rect_methy+rect_unmethy, r=args.rate)
#     return pval


# every call umr for the region, normalize the counts scales according to the fdr
def rescale_rect_counts(rect_methy, rect_unmethy):
    #J
    if rect_methy.quantile(0.25)==0:
        scale_factor = 1.0/50
    else:
        scale_factor = 1.0/rect_methy.quantile(0.25)
    #J
    #scale_factor = 1.0/rect_methy.quantile(0.25)
    rect_methy = scale_factor * rect_methy
    rect_unmethy = scale_factor * rect_unmethy
    return rect_methy, rect_unmethy

# recompute pval for the umr region using whole unmethy and methy reads counts
# as well as considering the umr width's
def calculate_region_pval(df_umr, args):
    rect_methy = df_umr.loc[:,'total_methy_counts'] / df_umr.loc[:,'region_width']
    rect_unmethy = df_umr.loc[:,'total_unmethy_counts'] / df_umr.loc[:,'region_width']
    rect_methy, rect_unmethy = rescale_rect_counts(rect_methy, rect_unmethy)
    pval = pbinom(rect_methy, rect_methy+rect_unmethy, r=args.rate)
    return pval

# call umrs and could handle the input samples with/without replicates
def call_umr_multi(chr_merge, ch, args):
    df = chr_merge.copy()
    ncols = df.shape[1]
    df.columns = ['methy','unmethy'] * int(ncols/2)
    if ncols > 2:
        df['methy_all'] = df['methy'].sum(axis=1)
        df['unmethy_all'] = df['unmethy'].sum(axis=1)
    else:
        df['methy_all'] = df['methy']
        df['unmethy_all'] = df['unmethy']
    df['pos'] = df.index.values
    # add the index
    df.index = range(0,len(df))
    # store the backup df before sliding
    df_backup = df.copy()
    # sliding the window
    df['methy_all'] = np.ceil(df['methy_all'].rolling(
            window=args.window_len,win_type='triang',center=True).sum()/args.window_len)
    df['unmethy_all'] = np.ceil(df['unmethy_all'].rolling(
            window=args.window_len,win_type='triang',center=True).sum()/args.window_len)
    df['total_all'] = df['methy_all'] + df['unmethy_all']
    df = df.dropna()
    rate = 1.0 * (df['methy_all'].sum()) / (df['total_all'].sum())
    args.rate = rate
    # change the rate
    # rate = 0.5
    # reformat the column names to have a fixed order for columns
    # robject_failur: if rpy2 failed to be imported, the significant slower version of numpy will be used
    if robject_failure:
        disp('Fail to import rpy2.robjects and the computation speed will be much slower!')
        p_val = np.apply_along_axis(pbinom_slow, 1, df.loc[:,'methy_all'].values, df.loc[:,'total_all'].values, r=rate)
    else:
        p_val = pbinom(df.loc[:,'methy_all'], df.loc[:,'total_all'], r=rate)
    df['p_val'] = p_val
    df_sub = df[df.p_val <= args.pval] # this threshold is just for pre-select potential region
    # df_sub.index: is the relative position for the methylation sites
    # combine all the UMR together
    df_umr = BB_HMM_compute_multi(df_sub, df_backup, args)
    df_umr['chr'] = ch
    df_umr = df_umr.reindex(columns=['chr', 'pos_start', 'pos_end', 'index_start', 'index_end', 'umr_fdr', \
        'num_methy_sites', 'region_width', 'mean_methy_ratio', 'variance_methy_ratio', 'total_methy_counts', 'total_unmethy_counts' ])

    return df_umr


# combine the results across chroms together to a big dataframe
def combine_list_dfs_to_write(sample_all, chs):
    df_chs = [sample_all.get(ch,'') for ch in chs if any(sample_all.get(ch,''))]
    if len(df_chs) > 1:
        df_chs_write = df_chs[0].append(df_chs[1:])
    else:
        df_chs_write = df_chs[0]
    return df_chs_write


# formate dataframe df_out to bed format
def umr2bed(df_out, which_fdr, ch=None):
    df_bed = pd.DataFrame({'chr': df_out.chr,
                       'chromStart': df_out.pos_start.astype(np.int32),
                       'chromEnd': df_out.pos_end.astype(np.int32),
                       'name': df_out[which_fdr].values,
                       'score': 0,
                       'strand': '+',
                       'thickStart': df_out.pos_start.astype(np.int32),
                       'thickEnd': df_out.pos_end.astype(np.int32),
                       'itemRgb': 0
                       })
    df_bed = df_bed.reindex(columns=['chr', 'chromStart', 'chromEnd', 'name',
                                 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'])
    return df_bed


def write2bed(id2umr_filtered, outdir, id2sample, suffix='umr'):
    for i, file_name in enumerate(id2umr_filtered):
        df_bed = umr2bed(id2umr_filtered[i])
        df_bed.to_csv(outdir+id2sample[i]+'.'+suffix+'.bed', sep='\t', header=['# chr', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'], index=False)


# convert the methylation file to the visualization wig file
def bed2wig(file_name, outdir, SPAN=1, HEADER=None, chs=['chr1','chr2']):
    data = pd.read_csv(file_name, header = HEADER, sep='\t')
    data = data.iloc[:, [0, 1, 2, 6, 7]]
    data.columns = ['chr', 'pos', 'strand', 'methy', 'total']
    data.index = data.pos
    # modification on 12/28/2017 for merge the strand
    data.pos = data.pos + 1  # wig pos = bed pos + 1
    file_prefix = file_name.split('/')[-1].split('.')[0]
    # if the file exist and delte such files
    for suffix in ['_methy.wig','_unmethy.wig','.wig']:
        if os.path.exists(outdir+file_prefix+suffix):
            os.remove(outdir+file_prefix+suffix)

    for ch in chs:
        print('processing on {}'.format(ch))
        header = "variableStep chrom={0} span={1}".format(ch, SPAN)
        data_sub = data[data['chr'] == ch]
        data_sub = merge_strand(data_sub)
        data_sub['total'] = data_sub['methy'] + data_sub['unmethy']
        data_sub['ratio'] = data_sub['methy']/data_sub['total']
        data_sub['pos'] = data_sub.index
        data_wig = data_sub.reindex(columns=['pos','methy','unmethy','total','ratio'])
        data_wig.iloc[:, [0,1]].to_csv(outdir+file_prefix+'_methy.wig', sep='\t', header=[header, ''], index=False, mode='a')
        data_wig.iloc[:, [0,2]].to_csv(outdir+file_prefix+'_unmethy.wig', sep='\t', header=[header, ''], index=False, mode='a')
        data_wig.iloc[:, [0,4]].to_csv(outdir+file_prefix+'.wig', sep='\t', header=[header, ''], index=False, mode='a')
