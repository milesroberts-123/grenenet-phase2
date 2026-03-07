import msprime
import pyslim

# Source - https://stackoverflow.com/a/16713986
# Posted by Pierre GM, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-06, License - CC BY-SA 4.0
import argparse

# Define the parser
parser = argparse.ArgumentParser(description='Short sample app')

# Declare an argument (`--algo`), saying that the 
# corresponding value should be stored in the `algo` 
# field, and using a default value if the argument 
# isn't given
parser.add_argument('--ID', action="store", dest='ID', default=0)
parser.add_argument('-N', action="store", dest='N', default=0)
parser.add_argument('-L', action="store", dest='L', default=0)
parser.add_argument('-R', action="store", dest='R', default=0)

# Now, parse the command line arguments and store the 
# values in the `args` variable
args = parser.parse_args()

# Individual arguments can be accessed as attributes...
print(args.ID)
print(args.N)

# Do msprime sim
ts = msprime.sim_ancestry(samples=int(args.N), population_size=args.N,
sequence_length=args.L, recombination_rate=args.R)
slim_ts = pyslim.annotate(ts, model_type="WF", tick=1)
slim_ts.dump("msprime_results/" + args.ID + ".trees")
