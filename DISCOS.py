import os
import re
import pandas as pd
import configparser
import multiprocessing
import subprocess
from pathlib import Path

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')
input_folder = config['Paths']['input_folder']
output_folder = config['Paths']['output_folder']
dcs_folder = config['Paths']['dcs_folder']

# Ensure output directory exists
Path(output_folder).mkdir(parents=True, exist_ok=True)
Path(dcs_folder).mkdir(parents=True, exist_ok=True)

def run_discos(input_file):
    """Run DISCOS on a given input file and save the output."""
    try:
        file_name = os.path.basename(input_file).split('.')[0]
        output_file = os.path.join(output_folder, f"{file_name}_output.txt")
        output_csv_file = os.path.join(dcs_folder, f"{file_name}_DCS.csv")
        # Assuming DISCOS can be invoked via command line
        result = subprocess.run(
            ['./DISCOS', input_file],  # Replace with the actual command to run DISCOS
            capture_output=True,
            text=True
        )
        print(result.stdout)
        # Save the output to the specified output file
        with open(output_file, 'w') as f:
            f.write(result.stdout)
        # Save the DCS output to the specified CSV file
            # Extract the line containing "DCS:"
        for line in result.stdout.splitlines():
            if line.startswith("DCS:"):
                # Get the gene list by splitting after "DCS: "
                gene_list = line.split("DCS: ")[1]
                pd.DataFrame(gene_list.split(",")).to_csv(output_csv_file, index=False, header=False)

        print(f"Processed: {file_name}, saved to {output_file}")
    except Exception as e:
        print(f"Error processing {input_file}: {e}")

def main():
    # Get a list of all files in the input folder
    input_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.csv')]

    # Use multiprocessing to run DISCOS on each file
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(run_discos, input_files)

    # read the .txt output files, extract the patients, genes, time and size into a statistics.csv file
    # under the same folder as the output files
    statistics_file = os.path.join(output_folder, 'statistics.csv')
    # Prepare a dictionary to store statistics for each cancer type
    statistics = {}
    # Extract data from each output .txt file
    for file_name in os.listdir(output_folder):
        if file_name.endswith('_output.txt'):
            cancer_type = file_name.split('_output.txt')[0]
            file_path = os.path.join(output_folder, file_name)

            with open(file_path, 'r') as file:
                content = file.read()

                # Extract Patients, Genes, Time, and Size of DCS using regex
                patients = int(re.search(r'Patients:\s+(\d+)', content).group(1))
                genes = int(re.search(r'Genes:\s+(\d+)', content).group(1))
                runtime = round(float(re.search(r'Time \(seconds\):\s+([0-9.]+)', content).group(1)), 2)
                dcs_size = int(re.search(r'Size of DCS:\s+(\d+)', content).group(1))

                # Store the extracted information
                statistics[cancer_type] = [patients, genes, runtime, dcs_size]

    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame.from_dict(statistics, orient='index',
                                columns=['Patients', 'Genes', 'Runtime (s)', 'Size of DCS'])
    df = df.transpose()  # Transpose so that the columns are cancer types
    # Save to CSV
    df.to_csv(statistics_file)
    # remove the .txt files
    for file_name in os.listdir(output_folder):
        if file_name.endswith('_output.txt'):
            os.remove(os.path.join(output_folder, file_name))


if __name__ == '__main__':
    main()