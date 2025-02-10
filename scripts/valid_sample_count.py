import yaml
import csv

def count_valid_samples(yaml_string):
    # Parse YAML content
    data = yaml.safe_load(yaml_string)
    
    # Dictionary to store species counts
    species_counts = {}
    
    # Process each species
    for species, info in data['species_mapping'].items():
        # Skip commented out species (they won't be in the parsed data)
        if 'samples' in info:
            # Count valid samples (not commented out)
            valid_samples = [sample for sample in info['samples'] 
                           if isinstance(sample, str) and not sample.startswith('#')]
            species_counts[species] = len(valid_samples)
    
    # Sort by species name
    sorted_counts = dict(sorted(species_counts.items()))
    
    # Write to CSV
    with open('data/valid_species_counts.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['species', 'valid_samples'])
        for species, count in sorted_counts.items():
            writer.writerow([species, count])
    
    return sorted_counts

# Example usage:
with open('config.yaml', 'r') as file:
    yaml_content = file.read()
counts = count_valid_samples(yaml_content)

