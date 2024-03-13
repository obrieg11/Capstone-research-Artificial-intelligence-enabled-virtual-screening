import sys
import argparse
import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import gc  # Garbage collection
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(levelname)s: %(message)s')

def read_molecules_in_chunks(sdf_file, chunk_size=5000):
    chunk = []
    mol_block = ''
    with open(sdf_file, 'r') as file:
        for line in file:
            mol_block += line
            if line.strip() == '$$$$':
                mol = Chem.MolFromMolBlock(mol_block)
                if mol is not None:
                    chunk.append((mol, mol_block))
                mol_block = ''

                if len(chunk) == chunk_size:
                    yield chunk
                    chunk = []
        if chunk:
            yield chunk

def calculate_similarity(query_fp, library_fp):
    return DataStructs.TanimotoSimilarity(query_fp, library_fp)

def get_compound_id(mol_block):
    lines = mol_block.split('\n')
    for i, line in enumerate(lines):
        if line.strip().startswith('>  <IDNUMBER>'):
            return lines[i + 1].strip()  # Return the line following the tag
    return "Unknown ID"

def generate_fingerprint(mol, fp_type):
    if fp_type == 'ECFP4':
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2)  # ECFP4 equivalent
    elif fp_type == 'FCFP4':
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True)  # FCFP4 equivalent
    else:
        raise ValueError("Unsupported fingerprint type")

def batch_insert(cursor, data):
    cursor.executemany('INSERT INTO results VALUES (?, ?, ?, ?)', data)

def main(query_sdf, library_sdf, output_txt, output_sdf, fp_type, chunk_size):
    conn = sqlite3.connect(':memory:')
    cursor = conn.cursor()
    cursor.execute('CREATE TABLE results (compound_id TEXT, similarity REAL, smiles TEXT, mol_block TEXT)')

    query_mol, _ = next(read_molecules_in_chunks(query_sdf, chunk_size=1))[0]
    query_fp = generate_fingerprint(query_mol, fp_type)
    del query_mol  # Free up memory

    mol_counter = 0

    for chunk in read_molecules_in_chunks(library_sdf, chunk_size):
        batch_data = []
        for lib_mol, lib_mol_block in chunk:
            try:
                lib_fp = generate_fingerprint(lib_mol, fp_type)
                sim = calculate_similarity(query_fp, lib_fp)
                compound_id = get_compound_id(lib_mol_block)
                smiles = Chem.MolToSmiles(lib_mol)
                batch_data.append((compound_id, sim, smiles, lib_mol_block))
            except Exception as e:
                logging.error(f"Error processing molecule: {e}")

        if batch_data:
            batch_insert(cursor, batch_data)
            mol_counter += len(batch_data)
            logging.info(f"Processed {mol_counter} molecules")

        gc.collect()  # Explicit garbage collection after processing each chunk

    cursor.execute('SELECT * FROM results ORDER BY similarity DESC')
    conn.commit()  # Committing the transaction after all data is inserted
    with open(output_txt, 'w') as out_txt, open(output_sdf, 'w') as out_sdf:
        for row in cursor:
            out_txt.write(f"{row[0]}: Similarity={row[1]}, SMILES={row[2]}\n")
            out_sdf.write(row[3] + '\n$$$$\n')
    logging.info("Completed processing")
    conn.close()

if __name__ == '__main__':    
    parser = argparse.ArgumentParser(description="Calculate molecular similarities.")
    parser.add_argument("query_sdf", help="Path to the query molecule SDF file")
    parser.add_argument("library_sdf", help="Path to the library SDF file")
    parser.add_argument("output_txt", help="Path to the output text file")
    parser.add_argument("output_sdf", help="Path to the output SDF file")
    parser.add_argument("--fptype", choices=['ECFP4', 'FCFP4'], default='ECFP4', help="Type of fingerprint to use (ECFP4 or FCFP4)")
    parser.add_argument("--chunksize", type=int, default=10000, help="Number of molecules to process at a time")
    args = parser.parse_args()
    main(args.query_sdf, args.library_sdf, args.output_txt, args.output_sdf, args.fptype, args.chunksize)

