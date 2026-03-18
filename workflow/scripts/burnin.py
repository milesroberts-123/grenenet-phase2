import msprime
import pyslim
import tskit
import numpy as np
# Source - https://stackoverflow.com/a/16713986
# Posted by Pierre GM, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-06, License - CC BY-SA 4.0
import argparse

# Define the parser
parser = argparse.ArgumentParser(description='recapitate ')

# Declare an argument (`--algo`), saying that the 
# corresponding value should be stored in the `algo` 
# field, and using a default value if the argument 
# isn't given
parser.add_argument('--ID', action="store", dest='ID', default=0)
parser.add_argument('-N', action="store", dest='N', default=0)
parser.add_argument('-L', action="store", dest='L', default=0)
parser.add_argument('-R', action="store", dest='R', default=0)
parser.add_argument('--mu', action="store", dest='mu', default=0)
parser.add_argument('--tau', action="store", dest='tau', default=0)

# Now, parse the command line arguments and store the 
# values in the `args` variable
args = parser.parse_args()

# Individual arguments can be accessed as attributes...
print(args.ID)
print(args.N)

# load ts from slim
print("Loading ts from slim...")
slim_ts = tskit.load("slim_fst_results/" + args.ID + ".trees")

print(f"The tree sequence has {slim_ts.num_trees} trees\n"
      f"on a genome of length {slim_ts.sequence_length},\n"
      f"{slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,\n"
      f"and {slim_ts.num_mutations} mutations.")

individual_times = slim_ts.individuals_time
for t in np.unique(individual_times):
    print(f"There are {np.sum(individual_times == t)} individuals from time {t}.")

# recapitate ts
print("Recapitating...")
rts = pyslim.recapitate(slim_ts,
             recombination_rate=args.R,
             ancestral_Ne=int(args.N))

# add mutations to ts
print("Adding mutations...")
next_id = pyslim.next_slim_mutation_id(rts)

ts = msprime.sim_mutations(
  rts,
  rate=args.mu,
  model=msprime.SLiMMutationModel(type=0, next_id=next_id),
  keep=True,
)

# calculate fst between historical and modern
print("Calculate fst...")
result = ts.Fst(sample_sets=[ts.samples(time=(1,args.tau)), ts.samples(time=(args.tau, max(individual_times)))])
print(result)

with open('msprime_results/' + args.ID + ".txt", 'w') as f:
  f.write(f"{result}\n")

#
nodes = ts.tables.nodes
is_sample = nodes.flags & tskit.NODE_IS_SAMPLE

# sample node IDs
sample_nodes = [u for u in ts.samples()]

# filter by time
ancient_samples = [u for u in sample_nodes if nodes.time[u] > int(args.tau)]
modern_samples  = [u for u in sample_nodes if (nodes.time[u] <= int(args.tau)) and (nodes.time[u] > 0)]

# simplify
ancient_ts = ts.simplify(samples=ancient_samples)
modern_ts  = ts.simplify(samples=modern_samples)

ancient_n = int(ancient_ts.num_samples / 2)
ancient_names = [f"his_{i}indv" for i in range(ancient_n)]

modern_n = int(modern_ts.num_samples / 2)
modern_names = [f"mod_{i}indv" for i in range(modern_n)]

# write to vcf
with open("msprime_results/ancient_" + args.ID + ".vcf", "w") as f:
    ancient_ts.write_vcf(f, individual_names = ancient_names)

with open("msprime_results/modern_" + args.ID + ".vcf", "w") as f:
    modern_ts.write_vcf(f, individual_names = modern_names)

print("Done!")
