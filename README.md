# Selection of high redshift quasars (z > 4.5)

This repository is for the 4MOST/CHANGES project, focusing on high-redshift quasar selection using SED fitting and color-based cuts. The workflow is modular, with each major step in a dedicated folder.

# The selection flowchart

```
                    ┌─────────────────┐
                    │   Input Data    │
                    │   ~420 mil.     │
                    └─────────┬───────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │   Removing      │
                    │   Initial       │
                    │ Contaminants    │
                    └─────────┬───────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │   SED fitting   │
                    │   ~ 13 mil.     │
                    └─────┬───┬───────┘
                          │   │
                          │   ▼
                          │ ┌─────────────────┐
                          │ │   Color Cut     │
                          │ └─────────┬───────┘
                          │           │
                          ▼           ▼
                    ┌─────────────────┐
                    │ Statistical     │
                    │ Cuts            │
                    │ ~ 10 mil.       │
                    └─────────┬───────┘
                              │
                              ▼
                    ╱─────────────────╲
                   ╱  Prioritization   ╲
                  ╱    ~ 930,000        ╲
                 ╱_____________________╲
                              │
                              ▼
                    ╲─────────────────╱
                     ╲ Initial Output ╱
                      ╲  ~ 24,500    ╱
                       ╲_____________╱
                              │
                              ▼
                        ◆─────────────────◆
                       ╱ Crossmatch with   ╲
                      ╱  DECaLS DR10.      ╲
                     ╱   Detection of      ╲
                    ╱    g_decals?         ╲
                   ◆─────────────────────────◆
                            │       │
                          yes│       │no
                             │       │
                             ▼       ▼
                   ┌─────────────────┐    ╭─────────────────╮
                   │ SNR cut of      │    │ Sample with     │
                   │ g_decals        │    │ DELVE DR2       │
                   │ ~ 23,000        │    │ photometry      │
                   └─────────┬───────┘    │ 1108 sources    │
                             │            ╰─────────┬───────╯
                             ▼                      │
                   ╭─────────────────╮              │
                   │ Sample with     │              │
                   │ DECaLS DR10     │              │
                   │ photometry      │              │
                   │ 5017 sources    │              │
                   ╰─────────┬───────╯              │
                             │                      │
                             ▼                      │
                    ╱─────────────────╲              │
                   ╱     Merging       ╲◄────────────┘
                  ╱___________________╲
                              │
                              ▼
                    ┌─────────────────┐
                    │  Final catalog  │
                    │  6125 sources   │
                    └─────────────────┘
```

## Process Flow Description

1. **Input Data** (~420 million sources) → Initial dataset
2. **Removing Initial Contaminants** → Data cleaning step
3. **SED fitting** (~13 million) → Spectral energy distribution analysis
4. **Color Cut** + **Statistical Cuts** (~10 million) → Quality filtering
5. **Prioritization** (~930,000) → Target selection
6. **Initial Output** (~24,500) → Preliminary results
7. **Crossmatch with DECaLS DR10** → Decision point for detection
   - **Yes**: SNR cut of g_decals (~23,000) → DECaLS DR10 photometry (5017 sources)
   - **No**: DELVE DR2 photometry (1108 sources)
8. **Merging** → Combine both branches
9. **Final catalog** (6125 sources) → Final output

## Folder Structure & Purposes

- **removing_initial_contaminants/**  
  Contains scripts and submodules for the first step: removing initial contaminants from the input data. This is the starting point of the workflow. It containts test results on each of the criteria in a separate folder and sed_cut.py which is a script of the initial cuts.

- **sed_script/**  
  Contains all SED fitting scripts and their dependencies. This is where the SED fitting step is performed after contaminants are removed and the magnitudes are converted to fluxes (mJy).  
  - `input_test/` and `output_test/` subfolders are used for input and output data, respectively, to keep the logic and file paths consistent.

  **Key scripts:**

  - `sed_calculation_command.py`:  
    This script performs SED fitting on input catalogs.  
    **How to run:**  
    ```
    python sed_calculation_command.py -i <input_file>
    ```
    where `<input_file>` is a file (e.g., `F23_fluxes_test.dat`) located in the `input_test/` subfolder.  
    **Test command:**  
    - Example: python sed_calculation_command.py -i F23_fluxes_test.dat (High-z qso samples from Fan+2023)
    
    **Outputs:**
    - The main output is a results file saved in the `output_test/` subfolder. The output file contains SED fitting results for each object in the input catalog, including best-fit parameters, chi-squared values, and statistical metrics for each template.


  - `sed_plot.py`:  
    This script generates plots from the SED fitting results and template files.  
    **How to run:**  
    ```
    python sed_plot.py
    ```
    The script will prompt for the object name, input catalog name, and list type (1 for QSO, 2 for BD, 0 for test).   

      **Test run:**  
   - Input the name of object J000009.99-041626.09 
   - Input the name of catalog F23_fluxes_test.dat
   - Input the type of list (1 for QSO, 2 for BD, 0 for test): 1

    **Outputs:**
    - The script generates plots visualizing the SED fitting results and template comparisons. These plots are saved as image files (e.g., PNG) in the `output_test/` subfolder.

- **color_cut/**  
  Contains scripts for applying color-based cuts to the data after SED fitting. This step further refines the candidate selection.

- **paper/**  
  Contains scripts, plots, and results used in publications and documentation of the project (Mkrtchyan et al. 2025 in prep.) Provided only with access.

- **query/**  
  Contains scripts for querying external catalogs (e.g., Gaia, DELVE, DECaLS) to crossmatch and supplement the main dataset.

- **output_catalogs/**  
  Stores the final and intermediate output catalogs generated at various steps of the workflow. Provided only with access.

- **tests/**  
  Contains test scripts and data for validating the pipeline and its components.

- **qso_bd_samples/**  
  Contains reference samples of known quasars and brown dwarfs for validation and comparison.Provided only with access.

---

## Notes

- The workflow is modular; you can run each step independently.
- The schematic above reflects the logical flow and data reduction at each stage.
- For more details on each script, see the comments and docstrings within the code.

