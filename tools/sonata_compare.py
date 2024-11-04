import libsonata
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    prog='sonata_compare',
    description='Compare two Sonata files, (optimistically) check if they contain the same data.')

parser.add_argument('files', nargs=2, help="Sonata HDF5 files to compare")
parser.add_argument('-n', '--name', type=str, help='Population name', default="All")
args = parser.parse_args()

pop1 = libsonata.EdgeStorage(args.files[0]).open_population(args.name)
pop2 = libsonata.EdgeStorage(args.files[1]).open_population(args.name)

print(f"Population 1 size: {pop1.size}")
print(f"Population 2 size: {pop2.size}")
assert pop1.size == pop2.size, "Population sizes should be equal in both files"

sel = pop1.select_all()

print("Checking source node populations...")
assert pop1.source == pop2.source, f"{pop1.source} != {pop2.source}"
print("Checking source node IDs equality...")
np.testing.assert_equal(pop1.source_nodes(sel), pop2.source_nodes(sel))
print("Checking target node populations...")
assert pop1.target == pop2.target, f"{pop1.target} != {pop2.target}"
print("Checking target node IDs equality...")
np.testing.assert_equal(pop1.target_nodes(sel), pop2.target_nodes(sel))

attrs1 = pop1.attribute_names
attrs2 = pop2.attribute_names

np.testing.assert_equal(attrs1, attrs2)

for attr in attrs1:
    print(f"Checking {attr} equality...")
    np.testing.assert_allclose(pop1.get_attribute(attr, sel), pop2.get_attribute(attr, sel))

print("All OK")
