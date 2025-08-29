import csv
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import re
import sys # Import sys to exit gracefully

# Define paths (consider making these configurable, e.g., via arguments or a config file)
BASE_DIR = r"C:\Users\aszyk\PycharmProjects\Mutation_modelling_complete"
PDB_DIR = os.path.join(BASE_DIR, "pdb")
REGION_DIR = os.path.join(BASE_DIR, "region_select")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
CSV_MUTANT_DIR = os.path.join(RESULTS_DIR, "CSV_mutant_files")
MODELLING_MUTANT_LISTS_DIR = os.path.join(BASE_DIR, "Modelling mutant lists")
MCSM_DIR = os.path.join(RESULTS_DIR, "mCSM")
MAESTRO_DIR = os.path.join(RESULTS_DIR, "MAESTRO output")
AVERAGE_DIR = os.path.join(RESULTS_DIR, "AVERAGE")
RESIDUE_PLOTS_DIR = os.path.join(RESULTS_DIR, "residue_plots_all")

# Create necessary directories
def create_directories():
    """Creates all necessary output directories."""
    dirs_to_create = [
        CSV_MUTANT_DIR,
        MODELLING_MUTANT_LISTS_DIR,
        MCSM_DIR,
        MAESTRO_DIR,
        AVERAGE_DIR,
        RESIDUE_PLOTS_DIR
    ]
    for directory in dirs_to_create:
        os.makedirs(directory, exist_ok=True)
    print("✅ All necessary directories created.")

def get_user_inputs():
    """Gets PDB ID, chain ID, and region file name from the user."""
    pdb_id = input(f"Please enter the PDB ID saved to {PDB_DIR}: ").strip()
    if not pdb_id:
        print("❌ PDB ID cannot be empty. Exiting.")
        sys.exit() # Exit if input is empty

    chain_id = input("Select chain ID: ").strip().upper()
    if not chain_id:
        print("❌ Chain ID cannot be empty. Exiting.")
        sys.exit() # Exit if input is empty

    region_name = input(f"Name of region selection .txt file saved in {REGION_DIR}: ").strip()
    if not region_name:
        print("❌ Region file name cannot be empty. Exiting.")
        sys.exit() # Exit if input is empty

    return pdb_id, chain_id, region_name

def parse_region_file(region_list_path):
    """Parses the region selection file."""
    region = set()
    if not os.path.exists(region_list_path):
        print(f"⚠️ Region file not found: {region_list_path}")
        return None
    else:
        with open(region_list_path) as f:
            for token in f:
                token = token.strip()
                if not token: continue
                for part in token.split(','):
                    if '-' in part:
                        try:
                            a, b = map(int, part.split('-'))
                            region.update(range(a, b + 1))
                        except ValueError:
                            print(f"⚠️ Invalid range in region file: {part}")
                    else:
                        try:
                            region.add(int(part))
                        except ValueError:
                            print(f"⚠️ Invalid position in region file: {part}")
    print(f"✅ Parsed {len(region)} mutation positions from region file.")
    return region

def parse_pdb_for_native_residues(pdb_path, chain):
    """Parses the PDB file to find native residues with CA atoms in the specified chain."""
    resname3to1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', # Corrected PRO back to P
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    native = {}
    if not os.path.exists(pdb_path):
        print(f"⚠️ PDB file not found: {pdb_path}")
        return None
    else:
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("ATOM"): continue
                resn = line[17:20].strip().upper()
                ch = (line[21] or ' ').strip() or 'A'
                try:
                    resi = int(line[22:26])
                except ValueError:
                    continue
                atom = line[12:16].strip()
                if ch != chain or atom != 'CA':
                    continue
                if resn in resname3to1:
                    native[resi] = resname3to1[resn]

    print(f"✅ Found {len(native)} residues with CA atoms in chain {chain}.")
    return native

def generate_mutation_list(native_residues, region, chain):
    """Generates a list of single-point mutations for a given PDB and region."""
    output_path = os.path.join(CSV_MUTANT_DIR, "CSV_mutant_files.csv")
    aa1 = set("ACDEFGHIKLMNPQRSTVWY")

    if not native_residues:
        print("❌ No native residues found. Cannot generate mutations.")
        return None
    if not region:
        print("❌ No region positions provided. Cannot generate mutations.")
        return None

    mutations = []
    for i in sorted(native_residues):
        if i not in region:
            continue
        wt = native_residues[i]
        for mut in sorted(aa1 - {wt}):
            mutations.append((chain, i, wt, mut))

    with open(output_path, "w", newline='') as f:
        w = csv.writer(f)
        w.writerow(["chain", "position", "wt", "mut"])
        for row in mutations:
            w.writerow(row)

    # Summary
    print("\n📊 Mutation Generation Summary")
    print(f"🧬 Mutation region positions considered: {len(region)}") # Clarified summary
    print(f"🔍 Native residues parsed: {len(native_residues)}")
    print(f"🧪 Mutations generated: {len(mutations)}")
    print(f"📁 Output saved to: {output_path}")
    return output_path

def csv_to_maestro(input_file):
    """Converts a CSV mutation list to MAESTRO format."""
    output_file = os.path.join(MODELLING_MUTANT_LISTS_DIR, "MAESTRO_mutant.txt")
    if not os.path.exists(input_file):
        print(f"⚠️ Input CSV file not found for MAESTRO conversion: {input_file}")
        return None

    # Read the CSV file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Skip the header
    lines = lines[1:]

    # Transform each line
    formatted_lines = []
    for line in lines:
        try:
            chain, position, wt, mut = line.strip().split(",")
            formatted_lines.append(f"{wt}{position}.{chain}{{{mut}}}")
        except ValueError:
            print(f"⚠️ Skipping invalid line in CSV: {line.strip()}")
            continue

    # Write to output file
    with open(output_file, "w") as f:
        f.write("\n".join(formatted_lines))
    print(f"📁 MAESTRO format file saved to: {output_file}") # Improved print message
    return output_file

def csv_to_mcsm(input_file):
    """Converts a CSV mutation list to mCSM format."""
    output_file = os.path.join(MODELLING_MUTANT_LISTS_DIR, "mCSM_mutant.txt")
    if not os.path.exists(input_file):
        print(f"⚠️ Input CSV file not found for mCSM conversion: {input_file}")
        return None

    # Read the CSV file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Skip the header
    lines = lines[1:]

    # Transform each line
    formatted_lines = []
    for line in lines:
        try:
            chain, position, wt, mut = line.strip().split(",")
            formatted_lines.append(f"{chain} {wt}{position}{mut}")
        except ValueError:
            print(f"⚠️ Skipping invalid line in CSV: {line.strip()}")
            continue

    # Write to output file
    with open(output_file, "w") as f:
        f.write("\n".join(formatted_lines))
    print(f"📁 mCSM format file saved to: {output_file}") # Improved print message
    return output_file

def transform_and_calculate_ddg():
    """Transforms mCSM and MAESTRO outputs and calculates average ΔΔG."""
    mcsm_input_path = os.path.join(MCSM_DIR, "mCSM_results.txt")
    maestro_input_path = os.path.join(MAESTRO_DIR, "MAESTRO_results.csv")
    transform_output_path = os.path.join(MCSM_DIR, "mCSM_results_transformed.csv")
    output_average_path = os.path.join(AVERAGE_DIR, "average.csv")

    if not os.path.exists(mcsm_input_path):
        print(f"⚠️ mCSM results file not found: {mcsm_input_path}. Skipping DDG calculation.")
        return None
    if not os.path.exists(maestro_input_path):
        print(f"⚠️ MAESTRO results file not found: {maestro_input_path}. Skipping DDG calculation.")
        return None

    try:
        # transform txt. output to CSV
        # Load the tab-delimited txt file
        df_mcsm = pd.read_csv(mcsm_input_path, sep="\t")

        # Normalize headers
        df_mcsm.columns = df_mcsm.columns.str.strip().str.upper()

        # Build mutation in file2 style: WILD_RES + RES_POS + "." + CHAIN + "{" + MUT_RES + "}"
        df_mcsm["MUTATION"] = df_mcsm["WILD_RES"] + df_mcsm["RES_POS"].astype(str) + "." + df_mcsm["CHAIN"] + "{" + df_mcsm["MUT_RES"] + "}"

        # Keep only mutation and ΔΔG
        df_mcsm_transformed = df_mcsm[["MUTATION", "PRED_DDG"]].rename(columns={"PRED_DDG": "DeltaG_tool1"})

        # Save as CSV
        df_mcsm_transformed.to_csv(transform_output_path, index=False)

        print("✅ mCSM file transformation done. Saved as file: " + transform_output_path)

        # calculate average ddG
        # Load transformed first file
        df1 = pd.read_csv(transform_output_path)

        # Load second CSV (semicolon-delimited)
        df2 = pd.read_csv(maestro_input_path, sep=";")

        # Strip whitespace from all column names
        df1.columns = df1.columns.str.strip()
        df2.columns = df2.columns.str.strip()

        # Create MUTATION column in df2 (use the substitution column directly)
        df2["MUTATION"] = df2["substitution"]

        # Rename ddG_pred to match df1
        df2 = df2.rename(columns={"ddG_pred": "DeltaG_tool2"})

        # Merge on MUTATION
        # Use 'outer' merge to keep all mutations from both files, fill missing with NaN
        merged = pd.merge(df1, df2[["MUTATION", "DeltaG_tool2"]], on="MUTATION", how="outer")


        # Compute average ΔΔG, handling potential NaNs from outer merge
        merged["DeltaG_avg"] = merged[["DeltaG_tool1", "DeltaG_tool2"]].mean(axis=1)

        # Save final CSV
        merged.to_csv(output_average_path, index=False)
        print("✅ The total values and averages were saved as " + output_average_path)
        return output_average_path

    except Exception as e:
        print(f"❌ Error during DDG calculation and merging: {e}")
        return None


def plot_all(average_ddg_path):
    """Generates individual plots for each residue's predicted ΔΔG values."""
    if not average_ddg_path or not os.path.exists(average_ddg_path):
        print("⚠️ Average ΔΔG file not found. Cannot generate all plots.")
        return None

    df = pd.read_csv(average_ddg_path)

    if df.empty:
         print("⚠️ Average ΔΔG DataFrame is empty. Cannot generate all plots.")
         return None

    # ---- Clean up MUTATION names ----
    # From format like "A34.A{C}" -> "A34C"
    # Use .loc[] for assignment to avoid SettingWithCopyWarning
    df.loc[:, "MUTATION_Clean"] = df["MUTATION"].str.replace(r"\..", "", regex=True)  # remove ".A" (chain)
    df.loc[:, "MUTATION_Clean"] = df["MUTATION_Clean"].str.replace(r"[{}]", "", regex=True)  # remove curly braces

    # Extract residue identifier (wildtype + position)
    df.loc[:, "Residue"] = df["MUTATION_Clean"].apply(lambda x: re.match(r"([A-Z]\d+)", x).group(0) if re.match(r"([A-Z]\d+)", x) else None)
    df = df.dropna(subset=["Residue"]) # Drop rows where residue could not be extracted

    df_long = df.melt(
        id_vars=["MUTATION", "MUTATION_Clean", "Residue"], # Include original MUTATION and CleanMutation
        value_vars=["DeltaG_tool1", "DeltaG_tool2"],
        var_name="Method",
        value_name="ΔΔG"
    )

    # Rename tools
    df_long["Method"] = df_long["Method"].replace({
        "DeltaG_tool1": "mCSM",
        "DeltaG_tool2": "MAESTRO"
    })

    # Set Seaborn style for pastel colors + grid
    sns.set_theme(style="whitegrid")
    palette = sns.color_palette("pastel", 2)

    # Loop over residues and make one figure each
    for residue, subset in df_long.groupby("Residue"):
        plt.figure(figsize=(8, 5))
        sns.barplot(
            data=subset,
            x="MUTATION_Clean", # Use cleaned mutation name for x-axis
            y="ΔΔG",
            hue="Method",
            palette=palette
        )
        plt.xticks(rotation=60, ha="right") # Added ha="right" for better label alignment
        plt.title(f"Predicted Stability Changes for {residue}", fontsize=14, weight="bold")
        plt.ylabel("ΔΔG (kcal/mol)")
        plt.xlabel("Mutation")
        plt.legend(title="Method")
        plt.tight_layout()

        # Save each plot separately
        plt.savefig(
            os.path.join(RESIDUE_PLOTS_DIR, f"{residue}_ddG.png"),
            dpi=300
        )
        plt.close()
    print(f"✅ Clean residue-specific plots saved in {RESIDUE_PLOTS_DIR} folder.")

def plot_ddg_consensus(df, save_path):
    """Generates a consensus plot based on agreement and ΔΔG thresholds."""
    if df.empty:
        print("⚠️ Input DataFrame for consensus plot is empty.")
        return pd.DataFrame()

    # Ensure necessary columns exist after potential outer merge in DDG calculation
    required_cols = ["MUTATION", "DeltaG_tool1", "DeltaG_tool2"]
    if not all(col in df.columns for col in required_cols):
        print(f"❌ Missing required columns for consensus plot. Needed: {required_cols}")
        return pd.DataFrame()

    # Extract residue identifier (wildtype + position)
    # Use .loc[] for assignment to avoid SettingWithCopyWarning
    df.loc[:, "Residue"] = df["MUTATION"].apply(lambda x: re.match(r"([A-Z]\d+)", x).group(0) if pd.notna(x) and re.match(r"([A-Z]\d+)", x) else None)
    df = df.dropna(subset=["Residue"]) # Drop rows where residue could not be extracted

    # Clean up mutation labels -> A34.A{C} → A34C
    df.loc[:, "CleanMutation"] = df["MUTATION"].apply(lambda x: re.sub(r"([A-Z]\d+)\..*?\{([A-Z])\}", r"\1\2", str(x))) # Handle potential NaN in MUTATION

    # Rename tools
    df = df.rename(columns={
        "DeltaG_tool1": "mCSM",
        "DeltaG_tool2": "MAESTRO"
    })

    # Agreement in sign - handle NaN values which means one tool didn't have a prediction
    df["Sign_mCSM"] = df["mCSM"].apply(lambda x: "+" if pd.notna(x) and x > 0 else ("-" if pd.notna(x) and x < 0 else "0"))
    df["Sign_MAESTRO"] = df["MAESTRO"].apply(lambda x: "+" if pd.notna(x) and x > 0 else ("-" if pd.notna(x) and x < 0 else "0"))


    # Keep only mutations where both tools agree AND both tools have a value (sign is not "0")
    df_agree = df[(df["Sign_mCSM"] == df["Sign_MAESTRO"]) & (df["Sign_mCSM"] != "0")].copy()


    # Further filter: at least one |ΔΔG| > 0.75
    df_agree = df_agree[(df_agree["mCSM"].abs() > 0.75) | (df_agree["MAESTRO"].abs() > 0.75)]

    if df_agree.empty:
        print("⚠️ No mutations satisfy the filtering criteria for consensus plot.")
        return df_agree

    # Melt for plotting
    df_long = df_agree.melt(
        id_vars=["CleanMutation", "Residue"],
        value_vars=["mCSM", "MAESTRO"],
        var_name="Method",
        value_name="ΔΔG"
    )

    # Plot
    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=df_long,
        x="CleanMutation",
        y="ΔΔG",
        hue="Method",
        palette=["#8ecae6", "#ffb5a7"],  # pastel blue/pink
        edgecolor="black"
    )
    plt.xticks(rotation=60, ha="right") # Added ha="right" for better label alignment
    plt.title("Consensus Predictions (same sign, |ΔΔG| > 0.75)", fontsize=14, weight="bold")
    plt.ylabel("ΔΔG (kcal/mol)")
    plt.xlabel("Mutation")
    plt.legend(title="Method")
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"✅ Plot saved: {save_path}")

    return df_agree

def write_rosetta_mut_file(df_agree, pdb_chain, output_path=None):
    """Writes a mutation list in Rosetta format."""
    def to_rosetta_line(mut_str, pdb_chain):
        match = re.match(r"([A-Z])(\d+)([A-Z])", mut_str)
        if match:
            wt, pos, mut = match.groups()
            return f"{pdb_chain} {pos} {wt} {mut}"
        return None

    if df_agree.empty:
        print("⚠️ Input DataFrame for Rosetta file is empty.")
        return None

    # Use 'CleanMutation' column which is already in the desired format
    mutation_lines = df_agree["CleanMutation"].apply(lambda x: to_rosetta_line(x, pdb_chain)).dropna().tolist()

    if not mutation_lines:
         print("⚠️ No valid mutations found to write to Rosetta file.")
         return None


    if output_path is None:
        output_path = os.path.join(RESULTS_DIR, "mutations_for_rosetta.txt")

    with open(output_path, "w") as f:
        for line in mutation_lines:
            f.write(line + "\n")

    print(f"✅ Mutation file saved: {output_path}")
    return output_path

# --- Main execution flow ---
if __name__ == "__main__": # Added if __name__ == "__main__": guard
    create_directories()
    pdb_id, chain, region_name = get_user_inputs()

    pdb_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    region_list_path = os.path.join(REGION_DIR, f"{region_name}.txt")

    # Step 1: Make mutant list
    csv_mutant_list_path = None
    region_data = parse_region_file(region_list_path)
    native_residues = parse_pdb_for_native_residues(pdb_path, chain)

    if region_data is not None and native_residues is not None:
         csv_mutant_list_path = generate_mutation_list(native_residues, region_data, chain)
    else:
        print("❌ Failed to parse region or PDB file. Cannot proceed with mutation generation.")


    if csv_mutant_list_path:
        # Step 2: Convert to MAESTRO format
        maestro_mutant_list_path = csv_to_maestro(csv_mutant_list_path)

        # Step 3: Convert to mCSM format
        mcsm_mutant_list_path = csv_to_mcsm(csv_mutant_list_path)

        print("\n--- Next Step: Run External Tools ---")
        print("Please download the generated files:")
        if maestro_mutant_list_path:
            print(f"- MAESTRO: {maestro_mutant_list_path}")
        if mcsm_mutant_list_path:
             print(f"- mCSM: {mcsm_mutant_list_path}")

        print("\nUpload these files to the respective webservers (MAESTRO and mCSM).")
        print("Download the results and save them as:")
        print(f"- MAESTRO results: {os.path.join(MAESTRO_DIR, 'MAESTRO_results.csv')}")
        print(f"- mCSM results: {os.path.join(MCSM_DIR, 'mCSM_results.txt')}")
        input("\nOnce the results are saved, please press enter to continue ...")
    else:
        print("❌ Mutation list generation failed. Skipping format conversion and subsequent steps.")

    # Assuming the user has manually run the external tools and saved the results
    # Step 4: Transform results and calculate average ddG
    average_ddg_path = transform_and_calculate_ddg()

    if average_ddg_path:
        # Load the average ddG dataframe once here for use in multiple plots
        try:
            df_average = pd.read_csv(average_ddg_path)

            # Step 5: Plot all mutations per residue
            plot_all(average_ddg_path)

            # Step 6: Plot consensus mutations
            consensus_plot_path = os.path.join(RESULTS_DIR, "Final_consensus_ddG.png") # Simplified filename
            df_agree = plot_ddg_consensus(df_average, consensus_plot_path)

            # Step 7: Write Rosetta mutation file
            if not df_agree.empty:
                rosetta_mut_file_path = write_rosetta_mut_file(df_agree, pdb_chain=chain)
            else:
                print("⚠️ No consensus mutations found. Skipping Rosetta file generation.")

        except Exception as e:
            print(f"❌ Error loading or processing average DDG file: {e}")


    else:
        print("⚠️ Average ΔΔG calculation failed or results not found. Skipping plotting and Rosetta file generation.")