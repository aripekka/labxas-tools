# labxas-tools
Tools to optimize laboratory-scale X-ray absorption experiments

### Prerequests 
The Python programs work either on Python 2 or 3. Apart from standard libraries, the cross-platform library xraylib is needed. For installation instructions see https://github.com/tschoonj/xraylib/wiki/Installation-instructions. The programs themselves need no installing.

### Usage
Currently labxas-tools consist of two programs:
* **XAS_sample_optimizer.py** calculates the optimal amount of sample and measurement times based on the chemical formula
* **XAS_measurement_time_estimator** estimates the measurement time needed to obtain a specific statistical accuracy based on the measured countrates

The programs can be run from terminal with 'python program_name.py' or 'python3 program_name.py', or by './program_name.py' if 'Allow execution' flag is set. The programs ask user for the relevant user input.

### More info
For technical details, see the pdf in docs/.
