#!/usr/bin/env python3
import yaml
import os
from pathlib import Path

def create_softlinks(config_file, base_path="../diversitynrecombination/data/recombinationrates/"):
    # Create output directory if it doesn't exist
    output_dir = Path("data/rec_rates")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read config file
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    # Process each species
    for species, info in config["species_mapping"].items():
        reference = info["reference"]
        
        # Get the last part of the species name for the source directory
        species_parts = species.split("_")
        source_dir = species_parts[-1].lower()
        
        # Construct source path
        source_path = Path(base_path) / reference / source_dir / f"nonpar/{source_dir}_100000_recomb.bed"
        
        # Construct destination path
        dest_path = output_dir / f"{species}_rec_rates.tsv"
        
        # Create softlink if source exists
        if source_path.exists():
            # Remove existing link if it exists
            if dest_path.exists() or dest_path.is_symlink():
                dest_path.unlink()
            
            # Create relative symlink
            try:
                dest_path.symlink_to(os.path.relpath(source_path, dest_path.parent))
                print(f"Created symlink: {dest_path} -> {source_path}")
            except Exception as e:
                print(f"Error creating symlink for {species}: {e}")
        else:
            print(f"Warning: Source file not found for {species}: {source_path}")

if __name__ == "__main__":
    create_softlinks("config.yaml")