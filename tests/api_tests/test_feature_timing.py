# encoding: utf-8
from __future__ import print_function

import os
nthreads = 1
os.environ["MKL_NUM_THREADS"] = "{}".format(int(nthreads))
os.environ["NUMEXPR_NUM_THREADS"] = "{}".format(int(nthreads))
os.environ["OMP_NUM_THREADS"] = "{}".format(int(nthreads))



import pickle
import time

import numpy as np

from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core




from memory_profiler import memory_usage
import gc










def benchmark_sim(method="lu", step=10, radius=1, wavelengths=[600], verbose=False):
    gc.collect()
    
    ## ---------- Setup structure
#    geometry = structures.rectwire(step, L=10,H=2,W=5, mesh='cube')
    geometry = structures.sphere(step, radius)
    material = materials.dummy(2.0)
    #material = materials.gold(interpolate_order=3)
    n1, n2 = 1.0, 1.0  # constant environment
    struct = structures.struct(step, geometry, material, n1,n2, 
                               structures.get_normalization('cube'))
    
    ## ---------- Setup incident field
    field_generator = fields.planewave        # planwave excitation
    kwargs = dict(theta = [0.0])              # several polarizations
    efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)
    
    ## ---------- Simulation initialization
    sim = core.simulation(struct, efield)
    
    
    ## --- benchmark `core.scatter`
    t0 = time.time()
    memusage = memory_usage(
            (core.scatter, (sim,), dict(method=method))
                            )
    t1 = time.time()
    
    Ndp = len(sim.struct.geometry)
    delta_t = round((t1-t0)*1000.)
    peak_mem = round(max(memusage),2)
    
    if verbose:
        print ("N dipoles   : ", Ndp)
        print ("peak memory : ", peak_mem, "mb")
        print ("total time  : ", delta_t, "ms")
        print ()
    
    return Ndp, delta_t, peak_mem




methods = ['lu', 'cg', 'dyson', 'scipyinv', 'numpyinv', 'pinv2', 'superlu']
methods = ['lu', 'cuda']#, 'cg', 'dyson', 'scipyinv', 'numpyinv', 'pinv2', 'superlu']
max_dipoles = [2200] * len(methods)
radii   = ((np.arange(3, 1500.5, 24))**(1/3.))[::2][:8]

ndipoles_radii = []
for radius_test in radii:
    Ndp = len(structures.sphere(1, radius_test))
    print ("R=", radius_test, "  Ndipols=", Ndp)
    ndipoles_radii.append(Ndp)

print ()
print ("total number of radii:", len(radii))


#%%
results = {}
#try:
#    results = pickle.load(open("benchmark_results_nthreads{}.pcl".format(nthreads), "rb"))
#except:
#    print "no results file found, start new test."

for imethod, method in enumerate(methods):
    if method in results:
        print ("        ---- method '{}' already evaluated, skipping...".format(method))
        continue
    results[method] = []
    print ()
    print (" ------------------ testing method '{}':".format(method))
    
    for i, r in enumerate(radii):
        if ndipoles_radii[i] > max_dipoles[imethod]:
            continue
        
        print ("    --> sphere R={:.3f}...".format(r), end='')
        Ndp, delta_t, peak_mem = benchmark_sim(method=method, radius=r, verbose=0)
        results[method].append([Ndp, delta_t, peak_mem])
        print ("{}: Npd={}, done in {}s, mem={:.1f}MB".format(i+1, Ndp, round(delta_t/1E3, 1), peak_mem))
    
    results[method] = np.array(results[method]).T
    pickle.dump(results, open("benchmark_results_nthreads{}.pcl".format(nthreads), "wb"))




#%%
#print results
#print results['lu'][0]
pickle.dump(results, open("benchmark_results_nthreads{}.pcl".format(nthreads), "wb"))



if 1:
    
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(14,4.5))
    plt.subplot(121)
    plt.plot(results['lu'][0], results['lu'][1]/1E3, color='b')
    plt.xlabel('meshpoints')
    plt.ylabel('time per wavelength (s)', color='b')
    
    
    
    #plt.twinx()
    plt.subplot(122)
    plt.plot(results['lu'][0], results['lu'][2], color='r')
    plt.xlabel('meshpoints')
    plt.ylabel('total memory (MB)', color='r')
    
    
    plt.show()
    
    
