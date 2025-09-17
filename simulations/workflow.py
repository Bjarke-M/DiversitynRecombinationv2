from gwf import Workflow, AnonymousTarget
from templates import *
import itertools
gwf = Workflow()

nes=[1000] #100,1000,10000
replicates=[i for i in range(10)]
parameter_space = itertools.product(nes,replicates)


for ne, replicate in parameter_space:
    jobind_bgs_slim = f'slim_bgs_ne_{ne}_replicate_{replicate}'
    gwf.target_from_template(jobind_bgs_slim,
                             bgs_only(
                                out = f'results/bgs_ne_{ne}_replicate_{replicate}.tree',
                                out_mid = f'results/bgs_20000_ne_{ne}_replicate_{replicate}.tree',
                                ne = f'{ne}'
                                )
)


s = [0.1,0.01,0.001]
hh_nes=[100,1000,10000] 
hh_parameter_space = itertools.product(hh_nes, s, replicates)

for ne, s, replicate in hh_parameter_space:
    jobind_hh_slim = f'slim_hh_ne_{ne}_s_{s}_replicate_{replicate}'
    gwf.target_from_template(jobind_hh_slim,
                             hh_only(
                                out = f'results/hh_full_ne_{ne}_s_{s}_replicate_{replicate}.tree',
                                out_mid = f'results/hh_20000_ne_{ne}_s_{s}_replicate_{replicate}.tree',
                                ne = f'{ne}',
                                s = f'{s}'
                                )
)