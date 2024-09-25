import numpy as np
import os
import re

def datafile_conversion(input_file):
    # Extract the initial mass from the filename
    def extract_mass_from_filename(filename):
        pattern = r'(\d+p\d+)g'
        match = re.search(pattern, filename)
        if match:
            mass_str = match.group(1).replace('p', '.')
            return float(mass_str)
        return None

    # Function to extract the date from the filename
    def extract_date_from_filename(filename):
        pattern = r'_(\d{8})\.txt$'
        match = re.search(pattern, filename)
        if match:
            return int(match.group(1))  # Convert the date string to an integer (YYYYMMDD)
        return None

    # Determine if the file is a gas ramp
    def determine_is_gas(filename):
        return 'gas' in filename.lower()

    # Read the residual mass from the file content
    def read_residual_mass_from_file(filename):
        residual_mass = None
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith("Residual Mass (g):"):
                        residual_mass = float(line.split(":")[1].strip())
                        break
        except FileNotFoundError:
            print(f"Error: File {filename} not found.")
        return residual_mass

    # Combine results from each function
    initial_mass = extract_mass_from_filename(input_file)
    is_gas = determine_is_gas(input_file)
    residual_mass = read_residual_mass_from_file(input_file)

    # ************************

    T = np.zeros(10000)
    p = np.zeros(10000)
    p_atm_psig = 14.7  # atmospheric pressure in psi
    dt = 30.  # recording freq, seconds
    psi_2_Pa = 6894.76  # Pascals/psi
    V = 0.001065  # m^3
    R = 8.314  # J/molK
    p_atm = 101325.  # atm

    # Heater profile correction coefficients
    a_value = 2.585525517584984e-07
    b_value = -0.00047141985292227707
    c_value = 1.1198651953917218

    # Temperature correction value is  a(T)^2 + b(T) = c
    # T_K is the raw data temp in Kelvin
    # T_bar is T_K * correction value

    # Extract the date from the filename
    file_date = extract_date_from_filename(input_file)
    cutoff_date = 20240921  # September 21st, 2024

    # Read input data file
    with open(input_file, 'r') as f:
        next(f)  # Skip the header line
        n = 0
        while True:
            try:
                T[n], p[n] = [float(x) for x in f.readline(15).split()]
                n += 1
                next(f)
            except ValueError:
                break

    # Calculate temperature corrected initial values
    T_K0 = T[0] + 273.15
    cor_coeff = a_value * (T_K0 ** 2) + (b_value * T_K0) + c_value
    T_0 = cor_coeff * (T_K0)
    p_0 = p_atm + (p[0] * 10 - p_atm_psig) * psi_2_Pa


    output_file = os.path.join(folder, filename.replace('.txt', '.csv'))
    
    # Choose the correct formula based on the date
    with open(output_file, 'w') as f:
        if is_gas:
            f.write('Time (h),T (C),p (psig),T (K),p (Pa),T_bar (K), pT ratio\n')
        else:
            f.write('Time (h),T (C),p (psig),T (K),p (Pa),T_bar (K), W (g/mol)\n')

        for i in range(0, n):
            time = (i * dt) / 3600.0  # Create time domain
            T_K = T[i] + 273.15  # Convert to Kelvin
            cor_coeff = a_value * (T_K ** 2) + (b_value * T_K) + c_value  # Temperature correction coefficient formula
            T_bar = cor_coeff * T_K  # Temperature correction

            # Apply different pressure formula based on the file date
            if file_date and file_date < cutoff_date:
                p_1 = p_atm + (p[i] * 10 - p_atm_psig) * psi_2_Pa
                p_psig = p[i] * 10 - p_atm_psig  # Corrected pressure
            else:
                p_1 = p_atm + (p[i] - p_atm_psig) * psi_2_Pa
                p_psig = p[i] - p_atm_psig  # Corrected pressure

            if is_gas:
                W_bar = max(0, (initial_mass - residual_mass) / ((V / R) * (p_1 / T_bar - p_0 / T_0)))
                f.write(f"{time:.4f},{T[i]:.1f},{p_psig:.1f},{T_K:.1f},{p_1:.1f},{T_bar:.1f},{W_bar:.1f}\n")
            else:
                ratio_bar = (p_1 / p_0) * (T_0 / T_bar)
                f.write(f"{time:.4f},{T[i]:.1f},{p_psig:.1f},{T_K:.1f},{p_1:.1f},{T_bar:.1f},{ratio_bar:.4f}\n")

# Run the conversion on all files in the folder
folder = input("Please enter folder path: ")

# Function to check for existing CSV files and convert if necessary
for filename in os.listdir(folder):
    if filename.endswith('.txt'):
        txt_file = os.path.join(folder, filename)
        csv_file = os.path.join(folder, filename.replace('.txt', '.csv'))
        # Check if the CSV file already exists
        if not os.path.exists(csv_file):
            print(f"Converting {filename} to CSV...")
            datafile_conversion(txt_file)
        else:
            print(f"CSV file for {filename} already exists. Skipping...")