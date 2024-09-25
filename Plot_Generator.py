import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pandas as pd


def datafile_conversion(filename):
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
    initial_mass = extract_mass_from_filename(filename)
    is_gas = determine_is_gas(filename)
    residual_mass = read_residual_mass_from_file(filename)

    return {
        "initial_mass": initial_mass,
        "is_gas": is_gas,
        "residual_mass": residual_mass
    }
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

    # Extract the date from the filename
    file_date = extract_date_from_filename(input_file)
    cutoff_date = 20240921  # September 21st, 2024

    # We corrected the scale factor for all experiments after sept 21 2024

    # Trim arrays to the actual size
    T = T[:n]
    p = p[:n]

    # Create arrays
    time_array = [i * dt / 3600 for i in range(n)]
    data_pt = []
    p_1_list = []
    molecular_weights = []

    # Calculating pT Ratio and/or Molecular Weight
    for i in range(1, n):
        T_K = T[i] + 273.15
        cor_coeff = a_value * (T_K ** 2) + (b_value * T_K) + c_value
        T_bar = cor_coeff * T_K
        # Apply different pressure formula based on the file date
        if file_date and file_date < cutoff_date:
            p_1 = p_atm + (p[i] * 10 - p_atm_psig) * psi_2_Pa
            p_psig = p[i] * 10 - p_atm_psig  # Corrected pressure
        else:
            p_1 = p_atm + (p[i] - p_atm_psig) * psi_2_Pa
            p_psig = p[i] - p_atm_psig  # Corrected pressure
        p_1_list.append(p_1)
        pT_ratio = (p_1 / p_0) * (T_0 / T_bar)
        molecular_weight = max(0, (initial_mass - residual_mass) / (
                    (V / R) * (p_1 / T_bar - p_0 / T_0))) if residual_mass is not None else 0
        molecular_weights.append(molecular_weight)
        data_pt.append(pT_ratio)

    # Apply a moving average to smooth the temperature data
    window_size = 20  # Adjust this window size as needed
    T_smoothed = np.convolve(T, np.ones(window_size) / window_size, mode='valid')
    derivatives_T = np.diff(T_smoothed) / dt
    slope_threshold = 0.025  # Adjust this threshold based on your data
    steady_state_T = np.abs(derivatives_T) < slope_threshold
    time_array_adjusted = time_array[:len(steady_state_T) + 1]
    T_adjusted = T[:len(steady_state_T) + 1]

    # Suppress steady state detection for 30 minutes after a large temperature increase
    temp_increase_threshold = 40.0  # Threshold for large temperature increase (in Celsius)
    suppress_duration = int(1800 / dt)  # Suppress duration (in seconds) divided by dt to get number of data points
    suppress_flag = np.zeros(len(T), dtype=bool)

    for i in range(1, len(T)):
        if i > suppress_duration and T[i] - T[i - 40] > temp_increase_threshold:
            steady_state_T[i:i + suppress_duration] = False

    # Adjust steady state detection based on suppression flag
    for i in range(len(steady_state_T)):
        if suppress_flag[i]:
            steady_state_T[i] = False

    # Extend steady state periods with a buffer on the front end
    buffer_size = 12  # Number of points to extend at the front end
    extended_steady_state_T = np.copy(steady_state_T)

    # Implement buffer for steady state periods
    for i in range(len(steady_state_T)):
        if steady_state_T[i]:
            # Remove buffer from front end
            if i >= buffer_size and not np.any(steady_state_T[i - buffer_size:i]):
                extended_steady_state_T[i] = True
            # # Add buffer to back end
            # end = min(len(steady_state_T), i + buffer_size)
            # extended_steady_state_T[i:end] = True

    # Collect temperature and pressure data during steady state periods
    steady_state_T_data = []
    steady_state_p_data = []
    steady_state_mw_data = []

    for i in range(len(extended_steady_state_T)):
        if extended_steady_state_T[i]:
            steady_state_T_data.append(T_adjusted[i])
            steady_state_p_data.append(p[i])
            steady_state_mw_data.append(molecular_weights[i])

    steady_state_T_data = np.array(steady_state_T_data)
    steady_state_p_data = np.array(steady_state_p_data)
    steady_state_mw_data = np.array(steady_state_mw_data)

    # Calculate pT ratio for steady state data
    steady_state_pt_ratio = []
    for T_ss, p_ss in zip(steady_state_T_data, steady_state_p_data):
        T_K_ss = T_ss + 273.15
        cor_coeff_ss = a_value * (T_K_ss ** 2) + (b_value * T_K_ss) + c_value
        T_bar_ss = cor_coeff_ss * T_K_ss
        p_psig_ss = p_ss * 10 - p_atm_psig
        p_1_ss = p_atm + p_psig_ss * psi_2_Pa
        pt_ratio_ss = (p_1_ss / p_0) * (T_0 / T_bar_ss)
        steady_state_pt_ratio.append(pt_ratio_ss)
    steady_state_pt_ratio = np.array(steady_state_pt_ratio)

    # Plotting all graphs
    if initial_mass != 0:
        fig, ((ax1, ax3), (ax2, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(14, 15))
    else:
        fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # Plot temperature vs. time
    ax1.plot(time_array, T, color='orange', label='Temperature')
    ax1.set_ylabel('Temperature (C)')
    ax1.set_xlabel('Time (hours)')
    ax1.set_title('Temperature vs. Time with Steady State Periods', fontweight='bold')

    # Create a secondary y-axis for pressure
    ax1b = ax1.twinx()
    ax1b.plot(time_array[1:], p_1_list, color='lightblue', label='Pressure')
    ax1b.set_ylabel('Pressure (Pa)')

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower right')

    for i in range(len(extended_steady_state_T)):
        if extended_steady_state_T[i]:
            ax1.axvspan(time_array_adjusted[i], time_array_adjusted[min(i + 1, len(time_array_adjusted) - 1)],
                        color='palegreen', alpha=0.3,
                        label='Steady State' if i == 0 else "")

    # Plot pT ratio vs. time
    ax2.plot(time_array[1:], data_pt, label='pT Ratio')
    for i in range(len(extended_steady_state_T)):
        if extended_steady_state_T[i]:
            ax2.axvspan(time_array_adjusted[i], time_array_adjusted[min(i + 1, len(time_array_adjusted) - 1)],
                        color='palegreen', alpha=0.3,
                        label='Steady State' if i == 0 else "")
    ax2.set_ylabel('pT Ratio (pT current / pT initial)')
    ax2.set_xlabel('Time (hours)')
    ax2.set_title('pT Ratio Analysis', fontweight='bold')
    ax2.legend(loc='upper right')

    # Plot steady state temperature data
    ax3.plot(steady_state_T_data, color='orange', label='Steady State Temperature')
    ax3.set_ylabel('Temperature (C)')
    ax3.set_title('Steady State Temperature and Pressure', fontweight='bold')
    ax3.set_xlabel('Steady State Time (minutes)')

    # Create a secondary y-axis for pressure data
    ax3b = ax3.twinx()
    ax3b.plot(steady_state_p_data, color='lightblue', label='Steady State Pressure')
    ax3b.set_ylabel('Pressure (Pa)')

    # Combine legends from both axes
    lines3, labels3 = ax3.get_legend_handles_labels()
    lines3b, labels3b = ax3b.get_legend_handles_labels()
    ax3.legend(lines3 + lines3b, labels3 + labels3b, loc='lower right')

    # Plot steady state pT ratio
    ax4.plot(steady_state_pt_ratio, label='pT Ratio at Steady States')
    ax4.set_ylabel('pT Ratio (pT current / pT initial)')
    ax4.set_xlabel('Steady State Time (minutes)')
    ax4.set_title('Steady State pT Ratio', fontweight='bold')
    ax4.legend(loc='lower right')

    # Plot Molecular Weight
    if initial_mass != 0:

        ax5.plot(molecular_weights, label='Molecular Weight')
        ax5.set_ylabel('Molecular Weight)')
        ax5.set_xlabel('Time (hours)')
        ax5.set_title('Molecular Weight Analysis', fontweight='bold')
        ax5.legend()

        # Fixes y-axis range for visibility
        if np.any(steady_state_mw_data > 1e3):  # Adjust the threshold as needed based on your data
            # Calculate percentiles only if outliers are present
            percentile_low = np.percentile(steady_state_mw_data, 3)
            percentile_high = np.percentile(steady_state_mw_data, 97)

            # Plot steady state Molecular Weights
            ax6.plot(steady_state_mw_data, label='Molecular Weight at Steady States')
            ax6.set_ylabel('Molecular Weight')
            ax6.set_xlabel('Steady State Time (minutes)')
            ax6.set_title('Molecular Weight at Steady States', fontweight='bold')
            ax6.legend()

            # Set y-axis limits for steady state Molecular Weight plot
            ax6.set_ylim(percentile_low, percentile_high)
        else:
            # Plot without setting percentiles
            ax6.plot(steady_state_mw_data, label='Molecular Weight at Steady States')
            ax6.set_ylabel('Molecular Weight')
            ax6.set_xlabel('Steady State Time (minutes)')
            ax6.set_title('Molecular Weight at Steady States', fontweight='bold')
            ax6.legend()

    # Save the figure as PNG with the input file name
    input_filename = input_file.split('.')[0]  # Remove the extension
    plt.savefig(f'{input_filename}.png', dpi=300)  # Save with .png extension and adjusted filename, dpi as needed

    fig.suptitle(f'{input_file} Steady State Analysis', fontsize=16)
    plt.tight_layout()
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make room for the subtitle
    plt.show()