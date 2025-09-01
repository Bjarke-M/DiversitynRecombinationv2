from gwf import Workflow, AnonymousTarget
from gwf.executors import Conda
gwf = Workflow()

def bgs_only(out,param):
    inputs = []
    outputs = [out]
    options = {"cores" : 1, 'memory': "18g", 'walltime': "1-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={param} -d "OUT='{out}'" bgs_only.slim  
    '''.format(param=param,out=out)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


def hh_only(out, ne, s):
    inputs = []
    outputs = []
    options = {"cores" : 1, 'memory': "18g", 'walltime': "1-00:00:00", 'account': "primatediversity"}
    executor = Conda('slim')
    spec = f'''
    slim -d N={ne} -d S={s} -d "OUT='{out}'" hh_only.slim  
    '''.format(ne=ne, s=s,out=out)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)