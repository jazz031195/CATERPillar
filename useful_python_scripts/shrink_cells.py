import pandas as pd
import numpy as np

# Define your file paths here
input_csv = '/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/Original_voxel.csv'      # Replace with your actual input file name
output_csv = '/home/localadmin/Documents/MCDS/Permeable_MCDS/output/incoherent_blood_flow/State1.csv' # The name of the new file to be generated

def convert_active_to_rest(input_path, output_path):
    # 1. Load the CSV into a pandas DataFrame
    df = pd.read_csv(input_path, sep=' ')
    
    # 2. Create masks to identify the specific rows for each cell type
    mask_blood_vessel = df['cell_type'] == 'blood_vessel'
    mask_axon = df['cell_type'] == 'axon'
    
    # 3. Apply the transformation to blood vessels (divide by 1.5)
    # This reverses the 50% radius expansion (e.g., bringing 6 um back down to 4 um)
    columns_to_modify = ['inner_radius', 'outer_radius']
    df.loc[mask_blood_vessel, columns_to_modify] /= 1.5
    
    # 4. Apply the transformation to axons (divide by sqrt(1.01))
    # Since Volume is proportional to radius squared, dividing the radius by sqrt(1.01)
    # perfectly scales the volume down by exactly 1%.
    df.loc[mask_axon, columns_to_modify] /= np.sqrt(1.01)
    
    # 5. Save the modified DataFrame to a new CSV file
    df.to_csv(output_path, sep=' ', index=False)
    print(f"Successfully modified the geometries and saved to: {output_path}")

# Run the function
if __name__ == "__main__":
    convert_active_to_rest(input_csv, output_csv)