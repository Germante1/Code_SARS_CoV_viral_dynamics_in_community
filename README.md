# Code_SARS_CoV_viral_dynamics_in_community

Code and models associated with the article:  
**“The analysis of a large dataset of PCR-test results in the community reveals the impact of age, variants, and vaccination status on SARS-CoV-2 viral dynamics.”**

## Repository content

- **`simuArticle.R`**  
  Script used to simulate the datasets for the analyses.

- **`Simulation_databases/`**  
  Contains one simulated dataset for each scenario (**obs_S1, obs_S2, obs_S3.csv**) generated using `simuArticle.R`, as well as the corresponding individual viral-dynamics parameters (**param_S1, param_S2, param_S3_1.csv**) for the simulated populations.  
  The folder also includes a dataset of **20,000 individuals** (10,000 infected and 10,000 always negative), where the number of tests and sampling times were directly sampled from the real-world database. This dataset (**community_PCR_tests_dataset.csv**) is used to simulate scenario 3.

- **Monolix model**  
  - Structural model: **`Structural_model_Monolix.txt`**  
  - Project file: **`Monolix.mlxtran`**

- **Stan models**
  - **Models for infected individuals only:**  
    - `Model_infected.stan`  
    - Execution script: `Init_run_Stan_Model_infected.R`
  - **Models for the full population:**  
    - `Model_full_population_infection_status.stan`  
    - Execution script: `Init_run_Stan_Model_full_population_infection_status.R`

## Software versions

All scripts were created using the following versions:

- **R 4.3.3**  
- **rstan 2.32**  
- **Monolix Suite 2023R1**

## Contact

**Maxime BEAULIEU**  
maxime.beaulieu@inserm.fr
