import csv
import math

import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, MACCSkeys, rdDepictor
from rdkit.ML.Cluster import Butina

###User inputs###
csv_file = input('Provide a file path to a .csv file with two columns, one called ID and one called SMILES: ').strip().strip('"')
image_name = input('Provide a name for the image: ')
mol_per_rows = input('No. of molecules per row? ')
SMILES_yn = input('Align by SMILES? y/n ')
if SMILES_yn == 'y':
	SMILES_INPUT = input('SMILES string to align by? ')
sdf_yn = input('Do you want to generate a .sdf file of the compounds as well? (y/n) ')

###Getting structures into python###
with open(csv_file, newline="") as file:
	structures_file = csv.DictReader(file)

	###Making sure names of columns are correct###
	for col in structures_file.fieldnames:
		if 'ID' in col.upper():
			id_column = col
			break

	for col in structures_file.fieldnames:
		if 'SMILES' in col.upper():
			SMILES_column = col
			break

	fp_list = []

	for row in structures_file:
		title = row[id_column]
		
		smiles = row[SMILES_column]

		###Generating molecules and MACCS keys (for clustering)###
		mol_gen = Chem.MolFromSmiles(smiles)

		if mol_gen != None:

			maccs_fp=MACCSkeys.GenMACCSKeys(mol_gen)

			fp_list.append([title, smiles, mol_gen, maccs_fp])

	###Creating distance matrix (for clustering)###
	distance_matrix = np.zeros(((len(fp_list)),(len(fp_list))))
	
	for i in range(0, (len(fp_list))):
		for j in range(0, (len(fp_list))):
			fp_1 = fp_list[i][3]
			fp_2 = fp_list[j][3]

			similarity = DataStructs.TanimotoSimilarity(fp_1, fp_2)

			distance = 1 - similarity
			distance_matrix[i, j] = round(distance,3)

	###Clustering, makes image look far neater###
	cut_off = 0.4
	clusters = Butina.ClusterData(distance_matrix, len(fp_list), cut_off, isDistData=True)
	clusters = sorted(clusters, key=len, reverse=True)
	no_clusters = len(clusters)

	structure_list = []

	###Aligning the compounds, for neatness###
	if SMILES_yn == 'y':
		ref = Chem.MolFromSmiles(SMILES_INPUT)
		rdDepictor.Compute2DCoords(ref)

		###Making a list of compounds and DDDs, compounds are aligned###
		for i in range(0, (no_clusters)):
			for compound in clusters[i]:
				DDD = fp_list[compound][0]
				DDD_mol_gen = fp_list[int(compound)][2]
				AllChem.GenerateDepictionMatching2DStructure(DDD_mol_gen, ref, acceptFailure=True)
				structure_list.append(DDD)
				structure_list.append(DDD_mol_gen)

	###if not aligning by SMILES###
	else:
		for i in range(0, (no_clusters)):
			for compound in clusters[i]:
				DDD = fp_list[compound][0]
				DDD_mol_gen = fp_list[int(compound)][2]
				structure_list.append(DDD)
				structure_list.append(DDD_mol_gen)


	###For generating an image###
	img = Draw.MolsToGridImage([structure_list[i] for i in range(1, len(structure_list), 2)], molsPerRow=int(mol_per_rows),  subImgSize=(200,200), legends=[structure_list[i] for i in range(0, len(structure_list), 2)])
	
	###generating dimensions of image###
	width = int(mol_per_rows)*2
	length = math.ceil(len(fp_list)/int(mol_per_rows))*2

	###generating plot of image###
	plt.figure(figsize=(width, length))
	plt.imshow(img)
	plt.title(image_name, fontsize = 12)
	plt.axis('off')
	plt.show()

	###For generating a .SDF file###
	if sdf_yn == "y":

		sdf_image_name = image_name + ".sd"

		writer = Chem.SDWriter(sdf_image_name)

		for i in range(1, len(structure_list), 2):
			writer.write(structure_list[i])

		writer.close()





