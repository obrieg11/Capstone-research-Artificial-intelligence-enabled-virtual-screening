import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
# Set up argument parser
parser = argparse.ArgumentParser(description='Generate specified number of conformers for molecules from an SDF file using RDKit.')
parser.add_argument('--input_sdf', required=True, help='Input SDF file')
parser.add_argument('--output_sdf', required=True, help='Output SDF file')
parser.add_argument('--num_conformers', type=int, required=True, help='Number of conformers to generate for each molecule')
args = parser.parse_args()
# Read molecules from SDF
supplier = Chem.SDMolSupplier(args.input_sdf)
molecules = [mol for mol in supplier if mol is not None]
# Generate conformers
for mol in molecules:
    AllChem.EmbedMultipleConfs(mol, numConfs=args.num_conformers)
# Write conformers to output SDF
with Chem.SDWriter(args.output_sdf) as writer:
    for mol in molecules:
        for confId in range(mol.GetNumConformers()):
            writer.write(mol, confId=confId)
print(f"Conformers generated and saved to {args.output_sdf}")

