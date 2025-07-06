
To use NECom, you need to have a solver that can handle mixed integer nonlinear programming (MINLP). Currently supported solvers are:

- BARON
- BONMIN (requires OPTI toolbox on Windows)
- GUROBI

Ensure one of these solvers is installed and properly configured in your MATLAB environment. For BONMIN on Windows, download and install the [OPTI Toolbox](https://inverseproblem.co.nz/OPTI/).

## Installation

1. Clone or download this repository from [https://github.com/Jingyi-Cai/NECom](https://github.com/Jingyi-Cai/NECom).
2. Add all folders and their subfolders to your MATLAB path:
   - Open MATLAB.
   - Go to `Home` > `Set Path`.
   - Add the root directory of NECom and all its subfolders.
   - Save the path.

## Usage

### Running the Demonstrative Example

1. Navigate to the `example_pd` folder.
2. Open `example.m`.
3. At the beginning of the script, set your preferred solver (e.g., 'baron', 'bonmin', 'gurobi').
4. Run the script to perform NECom simulations.

### Reproducing Results for Case 4 (Snowdrift Game)

1. Navigate to `./example_snowdrift/invertaseGame`.
2. Run `runExozyme.m` to reproduce the results for the snowdrift game case study.

## Data

All data related to the algae and yeast case study can be found in `DataFile_yeast_algae.xls`.

## References

For more information about the methodology and applications of NECom, please refer to the following paper:

- Cai, J., et al. (2020). NECom: Predicting Nash equilibria for microbial metabolic interactions. *Bioinformatics*, 36(24), 5649-5656. [https://academic.oup.com/bioinformatics/article/36/24/5649/6033579](https://academic.oup.com/bioinformatics/article/36/24/5649/6033579)
