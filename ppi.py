# This file is used to get the PPI connections.

from Bio.PDB import PDBParser, MMCIFParser, Selection, NeighborSearch
import os

def get_interface_residues(pdb_id, chain_ids, threshold):
    """
    Find residues of given chains that are parts of the interface residues within a given distance threshold.

    Args:
    pdb_id (str): The PDB ID of the structure.
    chain_ids (list): List containing two chain IDs.
    threshold (float): Maximum distance between residues to consider them as interacting.

    Returns:
    dict: Dictionary containing residue indices for each chain that are part of the interface.
          Returns empty lists if no interactions are found.
    """
    
    filename = ''
    extensions = ['.pdb', '.cif']
    for ext in extensions:
        if os.path.exists(f"{pdb_id.lower()}{ext}"):
            filename = f"{pdb_id.lower()}{ext}"
            break

    if not filename:
        raise FileNotFoundError(f"No .pdb or .cif file found for PDB ID '{pdb_id}'.")

    file_ext = filename.split('.')[-1]



    # Construct the file path
    filename = f"{pdb_id.lower()}.{file_ext}"

    # Select the appropriate parser
    if file_ext.lower() == 'pdb':
        parser = PDBParser(QUIET=False)
    elif file_ext.lower() == 'cif':
        parser = MMCIFParser(QUIET=False)
    else:
        raise ValueError("Unsupported file format. Use 'pdb' or 'mmcif'.")

    # Load the structure
    structure = parser.get_structure(pdb_id, filename)

    # print(vars(structure)) 
    '''
    {'level': 'S', '_id': '1XYZ', 'full_id': None, 'parent': None, 'child_list': [<Model id=0>], 
    'child_dict': {0: <Model id=0>}, 'xtra': {}, 
    'header': {'name': 'a common protein fold and similar active site in two distinct families of beta-glycanases', 
    'head': 'glycosyltransferase', 'idcode': '1XYZ', 'deposition_date': '1995-06-07', 'release_date': '1996-01-29', 
    'structure_method': 'x-ray diffraction', 'resolution': 1.4, 
    'structure_reference': ['h.souchon,s.spinelli,p.beguin,p.m.alzari crystallization and preliminary diffraction analysis catalytic domain of xylanase z from clostridium thermj.mol.biol. v. 235 1348 1994 issn 0022-2836 '], 
    'journal_reference': 'r.dominguez,h.souchon,s.spinelli,z.dauter,k.s.wilson,s.chauvaux,p.beguin,p.m.alzari a common protein fold and similar active site in two distinct families of beta-glycanases. nat.struct.biol. v. 2 569 1995 issn 1072-8368 7664125 10.1038/nsb0795-569 ', 
    'author': 'P.M.Alzari,S.Spinelli,R.Dominguez', 'compound': {'1': {'misc': '', 'molecule': '1,4-beta-d-xylan-xylanohydrolase', 'chain': 'a, b', 'synonym': 'endo-1\\,4-beta-xylanase z, xylanase xynz', 'ec_number': '3.2.1.8', 'ec': '3.2.1.8', 'engineered': 'yes'}}, 
    'source': {'1': {'misc': '', 'organism_scientific': 'clostridium thermocellum', 'organism_taxid': '1515', 'strain': 'ncib 10682', 'expression_system': 'escherichia coli', 'expression_system_taxid': '562', 'expression_system_plasmid': 'pct1214 (puc8)'}}, 'has_missing_residues': True, 
    'missing_residues': [{'model': None, 'res_name': 'PRO', 'chain': 'A', 'ssseq': 491, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'A', 'ssseq': 492, 'insertion': None}, {'model': None, 'res_name': 'GLN', 
    'chain': 'A', 'ssseq': 493, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'A', 'ssseq': 494, 'insertion': None}, {'model': None, 'res_name': 'ASP', 'chain': 'A', 'ssseq': 495, 'insertion': None}, {'model': None, 'res_name': 'VAL', 'chain': 'A', 'ssseq': 496, 'insertion': None}, {'model': None, 'res_name': 'GLN', 'chain': 'A', 'ssseq': 497, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'A', 'ssseq': 498, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'A', 'ssseq': 499, 'insertion': None}, {'model': None, 'res_name': 'ASN', 'chain': 'A', 'ssseq': 500, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'A', 'ssseq': 501, 'insertion': None}, {'model': None, 'res_name': 'SER', 'chain': 'A', 'ssseq': 502, 'insertion': None}, {'model': None, 'res_name': 'VAL', 'chain': 'A', 'ssseq': 503, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'A', 'ssseq': 504, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'A', 'ssseq': 505, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'A', 'ssseq': 506, 'insertion': None}, {'model': None, 'res_name': 'GLN', 'chain': 'A', 'ssseq': 507, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'A', 'ssseq': 508, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'A', 'ssseq': 509, 'insertion': None}, {'model': None, 'res_name': 'ILE', 'chain': 'A', 'ssseq': 510, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'A', 'ssseq': 511, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'A', 'ssseq': 512, 'insertion': None}, {'model': None, 'res_name': 'ILE', 'chain': 'A', 'ssseq': 513, 'insertion': None}, {'model': None, 'res_name': 'SER', 'chain': 'A', 'ssseq': 514, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'A', 'ssseq': 515, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'A', 'ssseq': 836, 'insertion': None}, {'model': None, 'res_name': 'TYR', 'chain': 'A', 'ssseq': 837, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'B', 'ssseq': 491, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'B', 'ssseq': 492, 'insertion': None}, {'model': None, 'res_name': 'GLN', 'chain': 'B', 'ssseq': 493, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'B', 'ssseq': 494, 'insertion': None}, {'model': None, 'res_name': 'ASP', 'chain': 'B', 'ssseq': 495, 'insertion': None}, {'model': None, 'res_name': 'VAL', 'chain': 'B', 'ssseq': 496, 'insertion': None}, {'model': None, 'res_name': 'GLN', 'chain': 'B', 'ssseq': 497, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'B', 'ssseq': 498, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'B', 'ssseq': 499, 'insertion': None}, {'model': None, 'res_name': 'ASN', 'chain': 'B', 'ssseq': 500, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'B', 'ssseq': 501, 'insertion': None}, {'model': None, 'res_name': 'SER', 'chain': 'B', 'ssseq': 502, 'insertion': None}, {'model': None, 'res_name': 'VAL', 'chain': 'B', 'ssseq': 503, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'B', 'ssseq': 504, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'B', 'ssseq': 505, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'B', 'ssseq': 506, 'insertion': None}, {'model': None, 'res_name': 'GLN', 'chain': 'B', 'ssseq': 507, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'B', 'ssseq': 508, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'B', 'ssseq': 509, 'insertion': None}, {'model': None, 'res_name': 'ILE', 'chain': 'B', 'ssseq': 510, 'insertion': None}, {'model': None, 'res_name': 'PRO', 'chain': 'B', 'ssseq': 511, 'insertion': None}, {'model': None, 'res_name': 'THR', 'chain': 'B', 'ssseq': 512, 'insertion': None}, {'model': None, 'res_name': 'ILE', 'chain': 'B', 'ssseq': 513, 'insertion': None}, {'model': None, 'res_name': 'SER', 'chain': 'B', 'ssseq': 514, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'B', 'ssseq': 515, 'insertion': None}, {'model': None, 'res_name': 'GLY', 'chain': 'B', 'ssseq': 836, 'insertion': None}, {'model': None, 'res_name': 'TYR', 'chain': 'B', 'ssseq': 837, 'insertion': None}], 'keywords': 'glycosyl hydrolase, xylanase, family f/10 of glycosyl hydrolases, clostridium thermocellum, glycosyltransferase', 'journal': 'AUTH   R.DOMINGUEZ,H.SOUCHON,S.SPINELLI,Z.DAUTER,K.S.WILSON,AUTH 2 S.CHAUVAUX,P.BEGUIN,P.M.ALZARITITL   A COMMON PROTEIN FOLD AND SIMILAR ACTIVE SITE IN TWOTITL 2 DISTINCT FAMILIES OF BETA-GLYCANASES.REF    NAT.STRUCT.BIOL.              V.   2   569 1995REFN                   ISSN 1072-8368PMID   7664125DOI    10.1038/NSB0795-569'}}

    '''
    
    # Extract the given chains from the protein structure.
    try:
        chain1 = structure[0][chain_ids[0]]
        chain2 = structure[0][chain_ids[1]]
    except KeyError:
        return {chain_ids[0]: [], chain_ids[1]: []} 

    # Find interacting residues
    atoms1 = Selection.unfold_entities(chain1, 'A')  # get all atoms from chain A at the atomic level
    atoms2 = Selection.unfold_entities(chain2, 'A')  # get all atoms from chain B at the atomic level
    ns = NeighborSearch(atoms1) #Return all atoms/residues/chains/models/structures that have at least one atom within radius of center. What entity level is returned (e.g. atoms or residues) is determined by level (A=atoms, R=residues, C=chains, M=models, S=structures).

    # Find pairs within the threshold
    interface_residues = {chain_ids[0]: set(), chain_ids[1]: set()}
    
    for atom in atoms2:
        # find the neighboring atoms at a residue level
        neighboring_residues = ns.search(atom.coord, 5.0, level='R')
        if neighboring_residues:
            for res in neighboring_residues:
                interface_residues[chain_ids[0]].add(res.get_full_id()[3][1]) # [3]: This accesses the fourth element of the tuple, which is the Residue ID tuple. [1]: This further accesses the second element of the Residue ID tuple, which is the sequence identifier of the residue.[1]: This further accesses the second element of the Residue ID tuple, which is the sequence identifier of the residue.
                interface_residues[chain_ids[1]].add(atom.get_parent().get_full_id()[3][1])
    
    
    
    # Convert sets to sorted lists in a more readable way
    interface_residues_dict = {}
    for chain_id, residues in interface_residues.items():
        if residues:
            sorted_residues = sorted(residues)
            interface_residues_dict[chain_id] = sorted_residues
        else:
            interface_residues_dict[chain_id] = []
    return interface_residues_dict
    
# Example usage:
residues = get_interface_residues("1XYZ", ["A", "B"], 5.0)
#residues2 = get_interface_residues("1XYZ", ["A", "B"], 5.0, "pdb")
print(residues)
#print(residues2)

# {'A': [25, 55, 62, 150, 153, 155, 156, 175, 178, 186, 189, 191, 206, 207, 281, 320, 393, 397, 676, 677, 680, 709, 710, 711, 712, 741, 744, 745, 746, 747, 749, 786, 787, 788, 790], 'B': [116, 217, 232, 252, 272, 282, 297, 319, 348, 351, 382, 391, 398, 400, 427, 535, 538, 539, 560, 561, 565, 566, 569, 570, 573, 578, 604, 605, 606, 607, 610, 806, 809]}

# {'A': [25, 55, 62, 150, 153, 155, 156, 175, 178, 186, 189, 191, 206, 207, 281, 320, 393, 397, 676, 677, 680, 709, 710, 711, 712, 741, 744, 745, 746, 747, 749, 786, 787, 788, 790], 'B': [116, 217, 232, 252, 272, 282, 297, 319, 348, 351, 382, 391, 398, 400, 427, 535, 538, 539, 560, 561, 565, 566, 569, 570, 573, 578, 604, 605, 606, 607, 610, 806, 809]}
