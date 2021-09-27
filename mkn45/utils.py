import numpy as np
import pandas as pd
import vireoSNP

def load_ground_truth(n_mut, af_pct, size=None, tree=None):
    if size is None and tree is None:
        ground = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct/clonal_var.csv')
    elif tree is not None:
        ground = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + tree + '/clonal_var.csv')
    elif size is not None:
        ground = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + size + '/clonal_var.csv')

    ground[['chrom', 'pos', 'ref', 'alt']] = ground['Unnamed: 0'].str.split('_', expand=True)
    ground['mgatk'] = ground['pos'] + ground['ref'] + '>' + ground['alt']
    
    return ground

def load_donor1_informative():
    ##load donor1_informative.csv
    mgatk_old = pd.read_csv('data/donor1_informative_variants.csv')
    mgatk_old[['ref', 'alt']] = mgatk_old.nucleotide.str.split('>', expand=True)
    mgatk_old['mquad_name'] = 'chrM_' + mgatk_old.position.astype(str) + '_' + mgatk_old.ref + '_' + mgatk_old.alt

    return mgatk_old

def process_mquad_data(n_mut, af_pct, size=None, tree=None):
    if size is None and tree is None:
        ground = load_ground_truth(n_mut, af_pct)
        mquad = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct/mquad/BIC_params.csv')
    elif tree is not None:
        ground = load_ground_truth(n_mut, af_pct, tree)
        mquad = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + tree + '/mquad/BIC_params.csv')
    elif size is not None:
        ground = load_ground_truth(n_mut, af_pct, size)
        mquad = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + size + '/mquad/BIC_params.csv')

    mquad['truth'] = mquad['0'].isin(ground['Unnamed: 0'])
    mquad.deltaBIC.replace(0, -100000, inplace=True)
    mgatk_old = load_donor1_informative()
    mquad['previously_identified'] = mquad['0'].isin(mgatk_old['Unnamed: 0'])
    mquad = mquad[mquad.previously_identified == False]

    return mquad

def process_mgatk_data(n_mut, af_pct, size=None, tree=None):
    if size is None and tree is None:
        ground = load_ground_truth(n_mut, af_pct)
        mgatk = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct/mgatk_in_mquad.csv')
    elif tree is not None:
        ground = load_ground_truth(n_mut, af_pct, tree)
        mgatk = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + tree + '/mgatk_in_mquad.csv')
    elif size is not None:
        ground = load_ground_truth(n_mut, af_pct, size)
        mgatk = pd.read_csv('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + size + '/mgatk_in_mquad.csv')

    mgatk_old = load_donor1_informative()
    
    mgatk['truth'] = mgatk['Unnamed: 0'].isin(ground['mgatk'])
    mgatk['previously_identified'] = mgatk['Unnamed: 0'].isin(mgatk_old['Unnamed: 0'])
    mgatk = mgatk[mgatk.previously_identified == False]
    mgatk = mgatk.fillna(0)
    
    return mgatk

def process_monovar_data(n_mut, af_pct, size=None, tree=None):
    if size is None and tree is None:
        ground = load_ground_truth(n_mut, af_pct)
        mpr = pd.read_table('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct/output_MPR.table')
    elif tree is not None:
        ground = load_ground_truth(n_mut, af_pct, tree)
        mpr = pd.read_table('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + tree + '/output_MPR.table')
    elif size is not None:
        ground = load_ground_truth(n_mut, af_pct, size)
        mpr = pd.read_table('data/' + str(n_mut) + '_clonal_variants_fixedAF' + str(af_pct) + 'pct_' + size + '/output_MPR.table')
        
    mpr['name'] = mpr['CHROM'] + '_' + mpr['POS'].astype(str) + '_' + mpr['REF'] + '_' + mpr['ALT']
    mpr.index = mpr.name
    mpr['truth'] = mpr.name.isin(ground['Unnamed: 0'])
    mquad = process_mquad_data(n_mut, af_pct)
    mquad.index = mquad['0']
    merged_with_mquad = mpr.merge(mquad['0'], left_index=True, right_index=True, how="outer").fillna(0)
    merged_with_mquad['truth'].replace(0, False, inplace=True)

    return merged_with_mquad

def load_vcf_data(vcf_path):
    ##Usage: Takes vcf as input, returns AD, DP, AF dataframes with index
    print("loading data from " + vcf_path + "...")

    cell_vcf = vireoSNP.load_VCF(vcf_path, biallelic_only=True)
    cell_dat = vireoSNP.vcf.read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])

    VAR = cell_vcf['variants']
    AD = pd.DataFrame(cell_dat['AD'].todense(), index=VAR)
    DP = pd.DataFrame(cell_dat['DP'].todense(), index=VAR)
    AF = (AD/DP).fillna(0)

    
    return AD, DP, AF