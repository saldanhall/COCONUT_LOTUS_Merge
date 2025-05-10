
## Merges LOTUS and COCONUT in a single and standardized plant-compound repository ###
####InchiKey RdKit
####Plant names with NCBI and WCVP


import pandas as pd

# -------- 1. Load the Data --------
# Load COCONUT dataset
coconut_db_long = pd.read_csv("/Users/saldanhaluizleonardo/Library/Mobile Documents/com~apple~CloudDocs/INDIGENOMICS/Garden and herbarium lists/COCONUT/coconut_db_long_with_ncbiid.csv")

# Load LOTUS dataset
lotus_db = pd.read_csv("/Users/saldanhaluizleonardo/Library/Mobile Documents/com~apple~CloudDocs/INDIGENOMICS/Garden and herbarium lists/230106_frozen_metadata.csv")

# -------- 2. Filter Only Plant Species --------
coconut_filtered = coconut_db_long[coconut_db_long['kingdom_ncbiid'] == 33090]

# LOTUS: Include both 'kingdom_taxid_ncbi == 33090' AND 'Archaeplastida'
lotus_filtered_kingdom = lotus_db[lotus_db['kingdom_taxid_ncbi'] == 33090]
lotus_filtered_archaeplastida = lotus_db[lotus_db['organism_taxonomy_02kingdom'] == "Archaeplastida"]

# Merge both filters to get the complete plant species list
lotus_filtered = pd.concat([lotus_filtered_kingdom, lotus_filtered_archaeplastida]).drop_duplicates()


print(lotus_filtered['organism_taxonomy_ncbiid'].head(10))


# -------- 3. Assign Prioritized Species IDs in LOTUS --------
lotus_filtered['ncbi_id_final_lotus'] = lotus_filtered.apply(
    lambda row: row['organism_taxonomy_ncbiid'] if pd.notna(row['organism_taxonomy_ncbiid']) else None, axis=1
)
lotus_filtered['wcvp_id_final_lotus'] = lotus_filtered.apply(
    lambda row: row['wcvp_accepted_id'] if pd.isna(row['organism_taxonomy_ncbiid']) and pd.notna(row['wcvp_accepted_id']) else None, axis=1
)


print(coconut_filtered['wcvp_accepted_id'].head(10))


print(lotus_filtered['organism_taxonomy_02kingdom'].unique())

lotus_filtered = lotus_filtered[
    ~lotus_filtered['organism_taxonomy_02kingdom'].fillna("").isin(['Metazoa', 'Fungi'])
]

print(lotus_filtered['organism_taxonomy_02kingdom'].unique())


# -------- 4. Assign Prioritized Species IDs in COCONUT --------
coconut_filtered['ncbi_id_final_coconut'] = coconut_filtered.apply(
    lambda row: row['species_ncbiid'] if pd.notna(row['species_ncbiid']) else None, axis=1
)
coconut_filtered['wcvp_id_final_coconut'] = coconut_filtered.apply(
    lambda row: row['wcvp_accepted_id'] if pd.isna(row['species_ncbiid']) and pd.notna(row['wcvp_accepted_id']) else None, axis=1
)

# -------- 5. Convert None to NaN (optional, for consistency) --------
lotus_filtered.loc[:, ['ncbi_id_final_lotus', 'wcvp_id_final_lotus']] = lotus_filtered[['ncbi_id_final_lotus', 'wcvp_id_final_lotus']].astype('float')
coconut_filtered.loc[:, ['ncbi_id_final_coconut', 'wcvp_id_final_coconut']] = coconut_filtered[['ncbi_id_final_coconut', 'wcvp_id_final_coconut']].astype('float')

# Print results
print("LOTUS: Assigned NCBI and WFO IDs")
print(lotus_filtered[['organism_name', 'ncbi_id_final_lotus', 'wcvp_id_final_lotus']].head(20))

print("COCONUT: Assigned NCBI and WFO IDs")
print(coconut_filtered[['organisms', 'ncbi_id_final_coconut', 'wcvp_id_final_coconut']].head())







# 1. Count unique IDs in 'ncbi_id_final_lotus' and 'wcvp_id_final_lotus'
unique_ncbi_lotus = lotus_filtered['ncbi_id_final_lotus'].nunique()
unique_wcvp_lotus = lotus_filtered['wcvp_id_final_lotus'].nunique()

print(f" Unique NCBI IDs in LOTUS: {unique_ncbi_lotus}")         #Unique NCBI IDs in LOTUS: 21546
print(f" Unique WCVP IDs in LOTUS: {unique_wcvp_lotus}")         #Unique WCVP IDs in LOTUS: 6349

# 2. Count unique IDs in 'ncbi_id_final_coconut' and 'wcvp_id_final_coconut'
unique_ncbi_coconut = coconut_filtered['ncbi_id_final_coconut'].nunique()
unique_wcvp_coconut = coconut_filtered['wcvp_id_final_coconut'].nunique()

print(f" Unique NCBI IDs in COCONUT: {unique_ncbi_coconut}")     #Unique NCBI IDs in COCONUT: 15422
print(f" Unique WCVP IDs in COCONUT: {unique_wcvp_coconut}")     # Unique WCVP IDs in COCONUT: 0

# 3. Count unique NCBI IDs across both datasets
unique_ncbi_combined = pd.concat([
    lotus_filtered['ncbi_id_final_lotus'],
    coconut_filtered['ncbi_id_final_coconut']
]).nunique()

print(f" Total unique NCBI IDs (LOTUS + COCONUT): {unique_ncbi_combined}")       #Total unique NCBI IDs (LOTUS + COCONUT): 24616


# 4. Count unique WCVP IDs across both datasets
unique_wcvp_combined = pd.concat([
    lotus_filtered['wcvp_id_final_lotus'],
    coconut_filtered['wcvp_id_final_coconut']
]).nunique()

print(f" Total unique WCVP IDs (LOTUS + COCONUT): {unique_wcvp_combined}")      #Total unique WCVP IDs (LOTUS + COCONUT): 6349

# 5. Sum unique NCBI and WCVP IDs
final_total_species = unique_ncbi_combined + unique_wcvp_combined
print(f" Final total unique species after harmonization: {final_total_species}")         #Final total unique species after harmonization: 30965

print(lotus_filtered.columns)
print(coconut_filtered.columns)

print(lotus_filtered['organism_taxonomy_02kingdom'].unique())



########Sparse SMILES

from rdkit import Chem
import pandas as pd

# Função para converter SMILES para InChIKey (14 primeiros caracteres)
def smiles_to_inchikey(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToInchiKey(mol)[:14]
    except:
        return None  # Se falhar, retorna None automaticamente

# Função para aplicar a conversão de SMILES para InChIKey em um DataFrame
def convert_smiles_to_inchikey(df, smiles_column, new_column):
    if smiles_column in df.columns:
        print(f"Convertendo SMILES de {smiles_column} para InChIKey...")
        df = df.dropna(subset=[smiles_column])  # Remover NaN antes
        df[new_column] = df[smiles_column].apply(smiles_to_inchikey)
        df = df.dropna(subset=[new_column])  # Remover valores nulos após a conversão
    else:
        print(f"Coluna '{smiles_column}' não encontrada no DataFrame.")
    return df

# Aplicar a conversão em COCONUT (canonical_smiles → InChIKey_short)
coconut_filtered = convert_smiles_to_inchikey(coconut_filtered, 'canonical_smiles', 'InChIKey_short')

# Aplicar a conversão em LOTUS (structure_smiles_2D → InChIKey_short)
lotus_filtered = convert_smiles_to_inchikey(lotus_filtered, 'structure_smiles_2D', 'InChIKey_short')

# Salvar os dados convertidos
coconut_filtered.to_csv("/Users/saldanhaluizleonardo/Library/Mobile Documents/com~apple~CloudDocs/INDIGENOMICS/Garden and herbarium lists/COCONUT/coconut_filtered.csv", index=False)
lotus_filtered.to_csv("/Users/saldanhaluizleonardo/Library/Mobile Documents/com~apple~CloudDocs/INDIGENOMICS/Garden and herbarium lists/lotus_filtered.csv", index=False)







print(lotus_filtered.columns)

print(coconut_filtered.columns)

##############################################################################################################
###################### INCHIKEY FIRST 14 CHARACTERS (after parsing smiles with rdkit) ######################
##############################################################################################################

# 1. Search ncbi_id_final_lotus in organism_taxonomy_ncbiid and retrieve structure_inchi
inchi_df_ncbi_lotus = lotus_filtered[['ncbi_id_final_lotus', 'InChIKey_short']].dropna()

# 2. Search wcvp_id_final_lotus in wcvp_accepted_id and retrieve structure_inchi
inchi_df_wfo_lotus = lotus_filtered[['wcvp_id_final_lotus', 'InChIKey_short']].dropna()

# 3. Search ncbi_id_final_coconut in species_ncbiid and retrieve standard_inchi
inchi_df_ncbi_coconut = coconut_filtered[['ncbi_id_final_coconut', 'InChIKey_short']].dropna()

# 4. Count unique standard_inchi in ncbi_id_final_coconut
unique_standard_inchi_ncbi_coconut = inchi_df_ncbi_coconut['InChIKey_short'].nunique()

# 5. Count unique structure_inchi in wcvp_id_final_lotus
unique_structure_inchi_wfo_lotus = inchi_df_wfo_lotus['InChIKey_short'].nunique()

# 6. Count unique structure_inchi in inchi_df_ncbi_lotus
unique_structure_inchi_ncbi_lotus = inchi_df_ncbi_lotus['InChIKey_short'].nunique()

# 7. Combine all unique InChI values and count distinct ones
all_unique_inchis = set(inchi_df_ncbi_coconut['InChIKey_short']) | \
                    set(inchi_df_ncbi_lotus['InChIKey_short']) | \
                    set(inchi_df_wfo_lotus['InChIKey_short'])

unique_inchi_combined = len(all_unique_inchis)

# Print results
print(f"✅ Unique standard_inchi in NCBI (Coconut): {unique_standard_inchi_ncbi_coconut}")
print(f"✅ Unique structure_inchi in WFO (LOTUS): {unique_structure_inchi_wfo_lotus}")
print(f"✅ Unique structure_inchi in NCBI (LOTUS): {unique_structure_inchi_ncbi_lotus}")
print(f"✅ Final total unique InChI values across all sources: {unique_inchi_combined}")   #Final total unique InChI values across all sources: 117042



print(coconut_filtered["standard_inchi_key"].head(30))
print(lotus_filtered["structure_inchikey"].head(30))  






##############################################################################################################
######################INCHIKEY FIRST 14 CHARACTERS (before parsing smiles with rdkit)######################
##############################################################################################################

# 1. Extract first 14 characters of structure_inchikey in LOTUS (NCBI)
inchi_df_ncbi_lotus = lotus_filtered[['ncbi_id_final_lotus', 'structure_inchikey']].dropna()
inchi_df_ncbi_lotus['structure_inchikey'] = inchi_df_ncbi_lotus['structure_inchikey'].str[:14]

# 2. Extract first 14 characters of structure_inchikey in LOTUS (WFO)
inchi_df_wfo_lotus = lotus_filtered[['wcvp_id_final_lotus', 'structure_inchikey']].dropna()
inchi_df_wfo_lotus['structure_inchikey'] = inchi_df_wfo_lotus['structure_inchikey'].str[:14]

# 3. Extract first 14 characters of standard_inchi_key in COCONUT (NCBI)
inchi_df_ncbi_coconut = coconut_filtered[['ncbi_id_final_coconut', 'standard_inchi_key']].dropna()
inchi_df_ncbi_coconut['standard_inchi_key'] = inchi_df_ncbi_coconut['standard_inchi_key'].str[:14]

# 4. Count unique InChIKeys (first 14 characters) in NCBI (Coconut)
unique_standard_inchi_ncbi_coconut = inchi_df_ncbi_coconut['standard_inchi_key'].nunique()

# 5. Count unique InChIKeys (first 14 characters) in WFO (LOTUS)
unique_structure_inchi_wfo_lotus = inchi_df_wfo_lotus['structure_inchikey'].nunique()

# 6. Count unique InChIKeys (first 14 characters) in NCBI (LOTUS)
unique_structure_inchi_ncbi_lotus = inchi_df_ncbi_lotus['structure_inchikey'].nunique()

# 7. Combine all unique InChIKeys and count distinct ones
all_unique_inchikeys = set(inchi_df_ncbi_coconut['standard_inchi_key']) | \
                       set(inchi_df_ncbi_lotus['structure_inchikey']) | \
                       set(inchi_df_wfo_lotus['structure_inchikey'])

unique_inchikey_combined = len(all_unique_inchikeys)

# Print results
print(f"✅ Unique standard_inchi_key in NCBI (Coconut): {unique_standard_inchi_ncbi_coconut}")
print(f"✅ Unique structure_inchikey in WFO (LOTUS): {unique_structure_inchi_wfo_lotus}")
print(f"✅ Unique structure_inchikey in NCBI (LOTUS): {unique_structure_inchi_ncbi_lotus}")
print(f"✅ Final total unique InChIKeys (first 14 characters) across all sources: {unique_inchikey_combined}")









# 1. Search ncbi_id_final_lotus in organism_taxonomy_ncbiid and retrieve structure_inchi
inchi_df_ncbi_lotus = lotus_filtered[['ncbi_id_final_lotus', 'structure_smiles_2D']].dropna()

# 2. Search wcvp_id_final_lotus in wcvp_accepted_id and retrieve structure_inchi
inchi_df_wfo_lotus = lotus_filtered[['wcvp_id_final_lotus', 'structure_smiles_2D']].dropna()

# 3. Search ncbi_id_final_coconut in species_ncbiid and retrieve standard_inchi
inchi_df_ncbi_coconut = coconut_filtered[['ncbi_id_final_coconut', 'canonical_smiles']].dropna()

# 4. Count unique standard_inchi in ncbi_id_final_coconut
unique_standard_inchi_ncbi_coconut = inchi_df_ncbi_coconut['canonical_smiles'].nunique()

# 5. Count unique structure_inchi in wcvp_id_final_lotus
unique_structure_inchi_wfo_lotus = inchi_df_wfo_lotus['structure_smiles_2D'].nunique()

# 6. Count unique structure_inchi in inchi_df_ncbi_lotus
unique_structure_inchi_ncbi_lotus = inchi_df_ncbi_lotus['structure_smiles_2D'].nunique()

# 7. Combine all unique InChI values and count distinct ones
all_unique_inchis = set(inchi_df_ncbi_coconut['canonical_smiles']) | \
                    set(inchi_df_ncbi_lotus['structure_smiles_2D']) | \
                    set(inchi_df_wfo_lotus['structure_smiles_2D'])

unique_inchi_combined = len(all_unique_inchis)

# Print results
print(f"✅ Unique standard_inchi in NCBI (Coconut): {unique_standard_inchi_ncbi_coconut}")
print(f"✅ Unique structure_inchi in WFO (LOTUS): {unique_structure_inchi_wfo_lotus}")
print(f"✅ Unique structure_inchi in NCBI (LOTUS): {unique_structure_inchi_ncbi_lotus}")
print(f"✅ Final total unique InChI values across all sources: {unique_inchi_combined}")



print(coconut_filtered["standard_inchi_key"].head(30))
print(lotus_filtered["structure_inchikey"].head(30))           #Final total unique InChIKeys (first 14 characters) across all sources: 117046


print(coconut_filtered.columns)