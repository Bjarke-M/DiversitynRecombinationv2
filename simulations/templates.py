from gwf import Workflow, AnonymousTarget
from gwf.executors import Conda
gwf = Workflow()

def bgs_only(out,out_mid,ne):
    inputs = []
    outputs = [out,out_mid]
    options = {"cores" : 1, 'memory': "512g", 'walltime': "7-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={ne} -d "OUT='{out}'" -d "OUT_MID='{out_mid}'" bgs_only.slim  
    '''.format(ne=ne,out=out)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


def hh_only(out, out_mid, ne, s):
    inputs = []
    outputs = [out,out_mid]
    options = {"cores" : 1, 'memory': "512g", 'walltime': "7-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={ne} -d S={s} -d "OUT='{out}'" -d "OUT_MID='{out_mid}'" hh_only.slim  
    '''.format(ne=ne, s=s,out=out)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)