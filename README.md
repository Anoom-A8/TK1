# ANUM N1: LU Factorization Strategies

The `anum_n1` folder contains Octave/MATLAB scripts for testing different LU factorization strategies on banded matrices. These scripts are intended to analyze the efficiency and accuracy of various LU strategies, including block LU, LU with pivoting, and recursive LU.

## Files

- **test_all_strategies.m**: Main script that runs all LU strategies on predefined matrix sizes and configurations.
- **test_block_lu_strategy.m**: Script to test the Block LU factorization strategy.
- **test_lu_pivot_strategy.m**: Script to test the LU factorization strategy with partial pivoting.
- **test_recur_fix.m**: Modified script to test recursive LU factorization with optimizations.
- **test_recur_lu_strategy.m**: Script to test the Recursive LU factorization strategy.

## Prerequisites

You will need Octave or MATLAB installed. The scripts are designed to be compatible with Octave, which is open-source and can be downloaded [here](https://www.gnu.org/software/octave/download.html).

## Running the Code

To run each script individually, you can use the Octave GUI or CLI:

1. **Using Octave GUI:**
   - Open Octave.
   - Navigate to the `anum_n1` folder.
   - Open and run any `.m` file, such as `test_all_strategies.m`, in the GUI.

2. **Using CLI:**
   - Ensure Octave is installed and added to your system path.
   - Open a terminal or command prompt.
   - Navigate to the `anum_n1` folder.
   - Run any script with the following command:
     ```bash
     octave filename.m
     ```
   - For example:
     ```bash
     octave test_all_strategies.m
     ```

## Output

The scripts will display various information, including:
- The matrix size (`N`), band width (`p`, `q`) used for the tests.
- Execution time of each LU factorization.
- Norm calculations to verify accuracy.
- Propagation error calculations.

## Reproducing Experimental Results

To obtain the same results as our experiments:
1. Run `test_all_strategies.m` or any other files you wish to experiment with, to execute all LU factorization strategies with predefined parameters.
2. Review the output for execution times, accuracy, and condition numbers.

For any issues or questions, please refer to our experimental notes or contact the team.

---

Happy experimenting!
