import dimod
from dwave.samplers import PlanarGraphSampler
from time import perf_counter

from dwave.samplers.random import sampler

def sq_lat(L, beta=1.0, num_reads=4):
    edgeset = {((i,j),(i,j+1)) for i in range(L) for j in range(L-1)} | {((i,j),(i+1,j)) for i in range(L-1) for j in range(L)}
    bqm = dimod.BinaryQuadraticModel(dimod.SPIN).from_ising({}, {e:-beta for e in edgeset})
    sampler = PlanarGraphSampler()
    start = perf_counter()
    sampleset = sampler.sample(bqm, num_reads=num_reads, beta=beta)
    end = perf_counter()
    print(f"Per sample time for {num_reads} samples L={L} beta={beta}: {(end - start)/num_reads}")
    start = perf_counter()
    sampleset = sampler.sample(bqm, num_reads=0, beta=beta, marginals=True)
    end = perf_counter()
    print(f"Per sample time for marginals L={L} beta={beta}: {(end - start)}")
    
if __name__=='__main__':
    for L in range(3,16):
        sq_lat(L, num_reads=10)

