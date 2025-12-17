from gwf import Workflow, AnonymousTarget
from gwf.executors import Conda
gwf = Workflow()


def neutral_only(out,out_mid,ne,mu,c):
    inputs = []
    outputs = [out,out_mid]
    options = {"cores" : 1, 'memory': "200g", 'walltime': "7-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={ne} -d "OUT='{out}'" -d mu={mu} -d C={c} -d "OUT_MID='{out_mid}'" slim_scripts/neutral_only.slim  
    '''.format(ne=ne,out=out)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


def bgs_only(out,out_mid,ne,mu,c):
    inputs = []
    outputs = [out,out_mid]
    options = {"cores" : 1, 'memory': "200g", 'walltime': "7-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={ne} -d "OUT='{out}'" -d mu={mu} -d C={c} -d "OUT_MID='{out_mid}'" slim_scripts/bgs_only.slim  
    '''.format(ne=ne,out=out)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


def hh_only(out, out_mid, ne, s, mu,c):
    inputs = []
    outputs = [out,out_mid]
    options = {"cores" : 1, 'memory': "200g", 'walltime': "7-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={ne} -d S={s} -d "OUT='{out}'" -d mu={mu} -d C={c} -d "OUT_MID='{out_mid}'" slim_scripts/hh_only.slim
    '''.format(ne=ne, s=s,out=out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# def recapitate_and_analyze(scenario, ne, s, n_replicates, recomb_ends, recomb_rates, out_csv):
#     # Build input file list based on scenario
#     if scenario == 'hh':
#         inputs = [f'results/hh_full_ne_{ne}_s_{s}_replicate_{i}.tree' for i in range(n_replicates)]
#     else:
#         inputs = [f'results/{scenario}_ne_{ne}_replicate_{i}.tree' for i in range(n_replicates)]

#     outputs = [out_csv]
#     options = {"cores": 4, 'memory': "64g", 'walltime': "12:00:00", 'account': "primatediversity"}
#     executor = Conda('recapitation')
#     spec = f'''
#     python recapitate_analyze.py {scenario} {ne} {s} {n_replicates} {recomb_ends} {recomb_rates} {out_csv}
#     '''
#     print(spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


def recapitate_and_analyze_new(scenario, ne, s, c, mu, mu_key, recomb_ends, recomb_rates, out_csv, input_list):
    inputs = [input_list]
    outputs = [out_csv]
    options = {"cores": 4, 'memory': "64g", 'walltime': "12:00:00", 'account': "primatediversity"}
    executor = Conda('recapitation')
    list_to_string = ' '.join(i for i in input_list)
    spec = f'''
    python recapitate_analyze_new.py {scenario} {ne} {s} {c} {mu} {mu_key} {recomb_ends} {recomb_rates} {out_csv} {list_to_string}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)