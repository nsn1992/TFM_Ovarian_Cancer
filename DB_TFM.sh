#!/bin/bash

# This script will execute Jupyter notebooks sequentially, generating excel files used as input in the next step. Some files used as input
# are already in the directory of execution. See readme in Github 

# Step 1: Execute the clinical notebooks. This step allow us to do a prefiltering of the initial clinical databases from the centers and hospitals.
echo "Running LP_clinical.ipynb..."
jupyter nbconvert --to notebook --execute LP_clinical.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "LP_clinical executed successfully"
else
    echo "Error executing LP_clinical"
    exit 1
fi

echo "Running OVE_clinical.ipynb..."
jupyter nbconvert --to notebook --execute OVE_clinical.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "OVE_clinical executed successfully"
else
    echo "Error executing OVE_clinical"
    exit 1
fi

echo "Running RVB_clinical.ipynb..."
jupyter nbconvert --to notebook --execute RVB_clinical.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "RVB_clinical executed successfully"
else
    echo "Error executing RVB_clinical"
    exit 1
fi

echo "Running MDA_clinical.ipynb..."
jupyter nbconvert --to notebook --execute MDA_clinical.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "MDA_clinical executed successfully"
else
    echo "Error executing MDA_clinical"
    exit 1
fi

# Step 2: Execute preunify notebooks. This allows unifying the notation of the different clinical variables into a common one for all centers or hospitals.
echo "Running LP_preunify.ipynb..."
jupyter nbconvert --to notebook --execute LP_preunify.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "LP_preunify executed successfully"
else
    echo "Error executing LP_preunify"
    exit 1
fi

echo "Running OVE_preunify.ipynb..."
jupyter nbconvert --to notebook --execute OVE_preunify.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "OVE_preunify executed successfully"
else
    echo "Error executing OVE_preunify"
    exit 1
fi

echo "Running MDA_preunify.ipynb..."
jupyter nbconvert --to notebook --execute MDA_preunify.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "MDA_preunify executed successfully"
else
    echo "Error executing MDA_preunify"
    exit 1
fi

echo "Running RVB_preunify.ipynb..."
jupyter nbconvert --to notebook --execute RVB_preunify.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "RVB_preunify executed successfully"
else
    echo "Error executing RVB_preunify"
    exit 1
fi

# Step 3: Execute the TILs_raw_counts notebook. This allows us to incorporate TILs raw counts information to the database containing the rest information
#  on genomic and genetic features and TILs scores (this database is also used as an input and it is called "Noelia Database complete.ods"
echo "Running TILs_raw_counts.ipynb..."
jupyter nbconvert --to notebook --execute TILs_raw_counts.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "TILs_raw_counts executed successfully"
else
    echo "Error executing TILs_raw_counts"
    exit 1
fi

# Step 4: Execute the final Combined_cohorts notebook. This script allow us to unify all type of variables (clinical,genomic,genetic and immunogenic)
# to generate a final database of work called "Samples_alltypedata_annotated.xlsx"
echo "Running Combined_cohorts_final.ipynb..."
jupyter nbconvert --to notebook --execute Combined_cohorts_final.ipynb --inplace
if [ $? -eq 0 ]; then
    echo "Combined_cohorts_final executed successfully"
else
    echo "Error executing Combined_cohorts_final"
    exit 1
fi

# Final message when all notebooks have been executed successfully
echo "All notebooks have been executed successfully."
