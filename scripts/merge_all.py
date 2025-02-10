import pandas as pd
import os

def extract_species_from_filename(filename):
    """
    Extracts species name from a path like 'data/stats/Alouatta_puruensis.combined.stats.csv'
    
    Args:
        filename: Path to the file
    
    Returns:
        Species name (e.g., 'Alouatta_puruensis')
    """
    # Get the base filename without path and extract species name
    base_name = os.path.basename(filename)
    species = base_name.split('.')[0]  # Takes everything before first period
    return species

def read_and_add_species(file_path):
    """
    Reads a CSV file and adds a species column based on the filename
    
    Args:
        file_path: Path to the CSV file
    
    Returns:
        DataFrame with added species column
    """
    # Read the file
    df = pd.read_csv(file_path)
    
    # Extract species name from filename
    species = extract_species_from_filename(file_path)
    
    # Add species column with the extracted name
    df['species'] = species
    
    return df

def concatenate_files_with_species(file_list):
    """
    Reads multiple files, adds species column to each, and concatenates them
    
    Args:
        file_list: List of file paths
    
    Returns:
        Concatenated DataFrame with species column
    """
    # Read each file and add species column
    dataframes = [read_and_add_species(file) for file in file_list]
    
    # Concatenate all dataframes
    concatenated_df = pd.concat(dataframes, ignore_index=True)
    
    return concatenated_df

# # Example usage:
# file_list = [
#     'data/stats/Alouatta_puruensis.combined.stats.csv',
#     'data/stats/Allochrocebus_lhoesti.combined.stats.csv'
# ]

# # Read files, add species column, and concatenate
# result = concatenate_files_with_species(file_list)

# # Save the result
# result.to_csv('test_concat_stats.csv', index=False)


file_list = snakemake.input.file_list
outfile = snakemake.output.all_stats

result = concatenate_files_with_species(file_list)
result.to_csv(outfile, index=False)