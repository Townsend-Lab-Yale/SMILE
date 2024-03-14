import csv
import numpy as np
import os
import sys
import shutil
from scipy.stats import entropy


def process_csv_file(csv_file_path, k_param=20):
    # Extract the filename from the path
    input_name = csv_file_path.split('/')[-1].split('.')[0]

    # Read CSV file
    with open(csv_file_path, 'r') as file:
        reader = csv.reader(file)
        
        # Calculate row and column numbers
        input = list(reader) # Convert reader to list
        row_number = len(input)
        file.seek(0)
        column_number = len(input[0])

        # Calculate stage and gene numbers
        stage_number = row_number - 2
        gene_number = column_number - row_number + 2


        # Assign parameter to k variable
        k = k_param

        # Gene_names
        header=next(reader) ## read the first row as the header
        gene_names = header[-gene_number:]


        # Print or use the calculated values as needed
        print(f"Step1: Read the input file")
        print(f"Row Number: {row_number}")
        print(f"Column Number: {column_number}")
        print(f"Stage Number: {stage_number}")
        print(f"Gene Number: {gene_number}")
        print(f"Parameter (k): {k}")

        # Return the reader and other calculated values
        return input_name, input, row_number, column_number, stage_number, gene_number, k, gene_names
    
def jensen_shannon_divergence(p, q):
    #Calculate Jensen-Shannon Divergence    
    m = (p + q) / 2  
    # the scipy entropy call is the Kullback-Leibler divergence. As scipy.stats.entropy normalizes the distributions, no need to normalize p and q
    return 0.5*(entropy(p,m)+entropy(q,m)) ## if not flatten 2-D array to 1-D array, the result is 1-D array, rather than a single value


def load_probability_matrix(file_path):
    #Load probability matrix from a file, skipping the header
    data = np.genfromtxt(file_path, skip_header=1, dtype=float, delimiter='\t', filling_values = np.nan) # filling_values=0, deletechars=''
    #return np.loadtxt(file_path, delimiter="\t", dtype=float, skiprows =1)
    # Remove empty rows and columns
    data = data[~np.isnan(data).all(axis=1)]
    data = data[:, ~np.isnan(data).all(axis=0)]
    return data

#def write_csv(output_matrix, output_file):
#    np.savetxt(output_file, output_matrix, delimiter=',', fmt='%f')

def write_csv(output_matrix, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerows(output_matrix)

if __name__ == '__main__':
     # Get file path and optional parameter (defaulting to 20) from command line arguments
    if len(sys.argv) < 2:
        print("Usage: python script.py <csv_file_path> [k]")
        sys.exit(1)
    
    csv_file_path = sys.argv[1]
    k_param = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    
    # Call the function to process the CSV file
    input_name, input_data, row_number, column_number, stage_number, gene_number, k, gene_names = process_csv_file(csv_file_path, k_param)
 
    #################################
    # Generate the run-k.sh file to specify the k value
    ##################################
    print(f"\nStep2: Generate the input directories and all necessary files for BNW")
    runscript = open("./tools/run.sh", "r")

    modi_run_file = open(f"./tmp/{input_name}_run-{k}.sh", "w")
    for line in runscript:
        if line != "k=100\n":
            modi_run_file.write(line)
        else:
            modi_run_file.write("k=" + str(k) + "\n")


    modi_run_file.close()
    runscript.seek(0) # Moves the file cursor to the beginning of the file
    
    ###############################
    # Generate the input dir and file
    ###############################
    input_data=np.matrix(input_data)
    #input_data = np.array(input_data, dtype = str)

    def generate_files(output_dir_prefix, output_dir_suffix, input_data = input_data,slice_columns=None):
        output_dir = f"{output_dir_prefix}{output_dir_suffix}"
        os.makedirs(output_dir, exist_ok=True)
        shutil.copy("./tools/banlist.txt", f"{output_dir}/banlist.txt")
        shutil.copy("./tools/whitelist.txt", f"{output_dir}/whitelist.txt")
        ## preserve the newlines within the file but avoid the extra newline at the end of the file.
        if slice_columns is not None:
            with open(f"{output_dir}/input.txt", "w",newline="") as file:
                for i, row in enumerate(input_data[:, slice_columns]):
                    file.write("\t".join(map(str, np.ravel(row))))
                    if i < len(input_data[:, slice_columns]) - 1:  # Add newline for all rows except the last one
                        file.write("\n")
           
        else:
            with open(f"{output_dir}/input.txt", "w", newline="") as file:
                for i, row in enumerate(input_data):
                    file.write("\t".join(map(str, np.ravel(row))))
                    if i < len(input_data) - 1:  # Add newline for all rows except the last one
                        file.write("\n")
          
    #for wild type with Gene and stage
    generate_files(f"./tmp/{input_name}_N0-0-{k}", "", input_data= input_data, slice_columns=None)
    #for wild type with GeneOnly
    generate_files(f"./tmp/{input_name}_N0-0-{k}","-GeneOnly", input_data= input_data,slice_columns=slice(stage_number, None))
   
    #for knock-out type
    for i in range(1, (gene_number + 1)): ## gene
        for j in range(1, (stage_number + 1)): ## stage
            cur = input_data.copy()
            cur[(j + 1):(stage_number + 2), (i + stage_number - 1)] = 0
            # for knock-out type with Gene and stage
            generate_files(f"./tmp/{input_name}_N{i}-{j}-{k}", "", input_data= cur, slice_columns=None)
            # for knock-out type with GeneOnly
            generate_files(f"./tmp/{input_name}_N{i}-{j}-{k}","-GeneOnly", input_data = cur, slice_columns=slice(stage_number, None))
    
    print(f"Generate the input directories and all necessary files for BNW: Finished!")


    #################################
    # run BNW to generate Bayesian Network, to get all the the model_averaging_probabilities.txt
    ##################################
    print(f"\nStep3: Run BNW to generate Bayesian Network")
     # run BNW for wild type with Gene and stage
    os.system(f"sh ./tmp/{input_name}_run-{k}.sh ./tmp/{input_name}_N0-0-{k} 1>./tmp/{input_name}_N0-0-{k}/run_BNW.log 2>&1")
    # run BNW for wild type with GeneOnly
    os.system(f"sh ./tmp/{input_name}_run-{k}.sh ./tmp/{input_name}_N0-0-{k}-GeneOnly 1>./tmp/{input_name}_N0-0-{k}-GeneOnly/run_BNW.log 2>&1")

    #for knock-out type
    for i in range(1, (gene_number + 1)): ## gene
        for j in range(1, (stage_number + 1)): ## stage
            # run BNW for knock-out type with Gene and stage
            os.system(f"sh ./tmp/{input_name}_run-{k}.sh ./tmp/{input_name}_N{i}-{j}-{k} 1>./tmp/{input_name}_N{i}-{j}-{k}/run_BNW.log 2>&1")
            # run BNW for knock-out type with GeneOnly
            os.system(f"sh ./tmp/{input_name}_run-{k}.sh ./tmp/{input_name}_N{i}-{j}-{k}-GeneOnly 1>./tmp/{input_name}_N{i}-{j}-{k}-GeneOnly/run_BNW.log 2>&1")

    print(f"Run BNW to generate Bayesian Network: Finished!")

    ##############################################
    # calculate the distance between Bayesian Networks, here we use Jensen Shannon Divergence
    ##############################################
    def calculate_jsd(input_name = input_name, gene_names = gene_names, gene_number=gene_number, stage_number=stage_number, k = k, input_type = ""):       
        # Wild type file path, and load probability matrix of wild type
        wild_type_file = f"./tmp/{input_name}_N0-0-{k}{input_type}/model_averaging_probabilities.txt"
        wild_type_matrix = load_probability_matrix(wild_type_file)
        #print(f"wild_type_matrix:{wild_type_file}")
        #print(wild_type_matrix)

        # Output matrix initialization
        output_matrix = np.zeros((stage_number, gene_number))
        rearranged_output_matrix = np.zeros((stage_number * gene_number, 2), dtype=object)

        # Loop through the knock-out type directories
        for gene in range(1, gene_number + 1):
            for stage in range(1, stage_number + 1):
                file_path = f"./tmp/{input_name}_N{gene}-{stage}-{k}{input_type}/model_averaging_probabilities.txt"

                # Load probability matrices
                current_matrix = load_probability_matrix(file_path)           
                #print(f"gene:{gene};stage:{stage}")
                #print(current_matrix)

                # Calculate Jensen-Shannon divergence
                jsd = jensen_shannon_divergence(wild_type_matrix.flatten(), current_matrix.flatten())

                # Store the result in the output matrix
                output_matrix[stage-1, gene-1] = jsd
                #print(f"gene:{gene};stage:{stage};jsd:{jsd}")

        # rearrange the matrix        
        for i in range(stage_number):
            for j in range(gene_number):            
                gene_name = f"{gene_names[j]}-{i + 1}"
                value = output_matrix[i, j]
                rearranged_output_matrix[i * gene_number + j, 0] = gene_name
                rearranged_output_matrix[i * gene_number + j, 1] = value

        # Print the output matrix 
        #print("Jensen-Shannon Divergence Matrix:")
        #print(output_matrix)
        #print("Rearranged Jensen-Shannon Divergence Matrix:")
        #print(rearranged_output_matrix)

        # Write the output matrix to a CSV file
        output_file = f"./output/{input_name}_jsd_matrix{input_type}.csv"
        output_matrix = np.vstack([gene_names, output_matrix])
        write_csv(output_matrix,output_file)
        print(f"Jensen-Shannon Divergence Matrix written to {output_file}")

        # Write the output matrix to a CSV file
        rearranged_output_file = f"./output/{input_name}_rearranged_jsd_matrix{input_type}.csv"
        write_csv(rearranged_output_matrix, rearranged_output_file)
        print(f"Rearranged Jensen-Shannon Divergence Matrix written to {rearranged_output_file}")
    
    print(f"\nStep4: Calculate the Jensen Shannon Divergence between Bayesian Networks")
    ## calculating JSD using the BNW outputs which input contains Gene and stage
    calculate_jsd(input_name = input_name, gene_names = gene_names, gene_number = gene_number, stage_number = stage_number, k= k, input_type="")
    
    ## calculating JSD using the BNW outputs which input contains Gene only
    calculate_jsd(input_name = input_name, gene_names = gene_names, gene_number = gene_number, stage_number = stage_number, k= k, input_type="-GeneOnly")

    print(f"\nDone!!!")