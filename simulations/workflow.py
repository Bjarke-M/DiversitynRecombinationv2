from gwf import Workflow, AnonymousTarget
from templates import *
import itertools
gwf = Workflow()

nes=[100,200,500,1000,1500,1800,2000]#,10000] 
c = [1,10,30,100]
s = [0.01]
mu = {'high':0.00000001}

replicates=[i for i in range(10)]
parameter_space = itertools.product(nes,c,mu.keys(),replicates)

for ne, c_val, mu_key, replicate in parameter_space:
    jobind_neutral_slim = f'slim_neutral_ne_{ne}_c_{c_val}_mu_{mu_key}_replicate_{replicate}'
    gwf.target_from_template(jobind_neutral_slim,
                             neutral_only(
                                out = f'results/neutral_ne_{ne}_mu_{mu_key}_c_{c_val}_replicate_{replicate}.tree',
                                out_mid = f'results/neutral_20000_ne_{ne}_mu_{mu_key}_c_{c_val}_replicate_{replicate}.tree',
                                ne = f'{ne}',
                                c = f'{c_val}',
                                mu = mu[f'{mu_key}']
                                ))
    jobind_bgs_slim = f'slim_bgs_ne_{ne}_c_{c_val}_mu_{mu_key}_replicate_{replicate}'
    gwf.target_from_template(jobind_bgs_slim,
                             bgs_only(
                                out = f'results/bgs_ne_{ne}_mu_{mu_key}_c_{c_val}_replicate_{replicate}.tree',
                                out_mid = f'results/bgs_20000_ne_{ne}_mu_{mu_key}_c_{c_val}_replicate_{replicate}.tree',
                                ne = f'{ne}',
                                c = f'{c_val}',
                                mu = mu[f'{mu_key}']
                                ))


parameter_space_hh = itertools.product(nes,c,mu.keys(),s,replicates)

for  ne, c_val, mu_key, s_val, replicate in parameter_space_hh:
    jobind_hh_slim = f'slim_hh_ne_{ne}_c_{c_val}_mu_{mu_key}_s_{s_val}_replicate_{replicate}'
    gwf.target_from_template(jobind_hh_slim,
                             hh_only(
                                out = f'results/hh_full_ne_{ne}_mu_{mu_key}_c_{c_val}_s_{s_val}_replicate_{replicate}.tree',
                                out_mid = f'results/hh_20000_ne_{ne}_mu_{mu_key}_c_{c_val}_s_{s_val}_replicate_{replicate}.tree',
                                ne = f'{ne}',
                                s = f'{s_val}',
                                c = f'{c_val}',
                                mu = mu[f'{mu_key}']
                                )
)

### Analysis targets - recapitate and analyze diversity-recombination relationship ###
analyze_space_hh = itertools.product(nes,c,mu.keys(),s)
for  ne, c_val, mu_key, s_val in analyze_space_hh:
# HH
    jobind_analyze_slim = f'analyze_hh_ne_{ne}_c_{c_val}_mu_{mu_key}_s_{s_val}'
    gwf.target_from_template(jobind_analyze_slim,
                             recapitate_and_analyze_new(
                                 scenario='hh', 
                                 ne=ne, 
                                 s=s_val, 
                                 c=c_val, 
                                 mu=mu[mu_key], 
                                 mu_key=mu_key, 
                                 recomb_ends = 'stdpopsim_chr12_ends1.csv',
                                 recomb_rates = 'stdpopsim_chr12_rates1.csv',
                                 out_csv = f'results/analyzed_hh_ne_{ne}_s_{s_val}_mu_{mu_key}_c_{c_val}.csv',
                                 input_list = [f'results/hh_full_ne_{ne}_mu_{mu_key}_c_{c_val}_s_{s_val}_replicate_{i}.tree' for i in range(10)]
                             ))
analyze_space = itertools.product(nes,c,mu.keys())
for  ne, c_val, mu_key in analyze_space:
# BGS
    jobind_analyze_slim = f'analyze_bgs_ne_{ne}_c_{c_val}_mu_{mu_key}'
    gwf.target_from_template(jobind_analyze_slim,
                             recapitate_and_analyze_new(
                                 scenario='bgs', 
                                 ne=ne, 
                                 s = 'NA',
                                 c=c_val, 
                                 mu=mu[mu_key], 
                                 mu_key=mu_key, 
                                 recomb_ends = 'stdpopsim_chr12_ends1.csv', 
                                 recomb_rates = 'stdpopsim_chr12_rates1.csv', 
                                 out_csv = f'results/analyzed_bgs_ne_{ne}_mu_{mu_key}_c_{c_val}.csv', 
                                 input_list = [f'results/bgs_ne_{ne}_mu_{mu_key}_c_{c_val}_replicate_{i}.tree' for i in range(10)]
                             ))
# Neutral 
    jobind_analyze_slim = f'analyze_neutral_ne_{ne}_c_{c_val}_mu_{mu_key}'
    gwf.target_from_template(jobind_analyze_slim,
                             recapitate_and_analyze_new(
                                 scenario='neutral', 
                                 ne=ne, 
                                 s = 'NA',
                                 c=c_val, 
                                 mu=mu[mu_key], 
                                 mu_key=mu_key, 
                                 recomb_ends = 'stdpopsim_chr12_ends1.csv', 
                                 recomb_rates = 'stdpopsim_chr12_rates1.csv', 
                                 out_csv = f'results/analyzed_neutral_ne_{ne}_mu_{mu_key}_c_{c_val}.csv', 
                                 input_list = [f'results/neutral_ne_{ne}_mu_{mu_key}_c_{c_val}_replicate_{i}.tree' for i in range(10)]
                             ))








# nes=[100,1000,10000]

# # Neutral simulations
# for ne in nes:
#     gwf.target_from_template(f'analyze_neutral_ne_{ne}',
#                              recapitate_and_analyze(
#                                 scenario='neutral',
#                                 ne=ne,
#                                 s='NA',
#                                 n_replicates=10,
#                                 recomb_ends='stdpopsim_chr12_ends1.csv',
#                                 recomb_rates='stdpopsim_chr12_rates1.csv',
#                                 out_csv=f'results/analyzed_neutral_ne_{ne}.csv'
#                                 )
# )

# # BGS simulations
# for ne in nes:
#     gwf.target_from_template(f'analyze_bgs_ne_{ne}',
#                              recapitate_and_analyze(
#                                 scenario='bgs',
#                                 ne=ne,
#                                 s='NA',
#                                 n_replicates=10,
#                                 recomb_ends='stdpopsim_chr12_ends1.csv',
#                                 recomb_rates='stdpopsim_chr12_rates1.csv',
#                                 out_csv=f'results/analyzed_bgs_ne_{ne}.csv'
#                                 )
# )

# # HH simulations
# s_values = [0.1, 0.01, 0.001]
# for ne in nes:
#     for s in s_values:
#         gwf.target_from_template(f'analyze_hh_ne_{ne}_s_{s}',
#                                  recapitate_and_analyze(
#                                     scenario='hh',
#                                     ne=ne,
#                                     s=s,
#                                     n_replicates=10,
#                                     recomb_ends='stdpopsim_chr12_ends1.csv',
#                                     recomb_rates='stdpopsim_chr12_rates1.csv',
#                                     out_csv=f'results/analyzed_hh_ne_{ne}_s_{s}.csv'
#                                     )
#         )
