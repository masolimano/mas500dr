# MAS500 Data Reduction Pipeline
Python-based utility to reduce imaging data taken with the CCD camera mounted on the MAS500 telescope on El Sauce Observatory.

## Installation (for development)
1. Create a Python 3.8 virtual environment for this project
```sh
python -m venv ~/.mfhenv --prompt mas500
source ~/.mfhenv/bin/activate
pip install setuptools wheel
```
2. Clone the project
```sh
git clone https://github.com/masolimano/mas500dr.git
cd mas500dr
```
3. Use pip to install it in development mode
```sh
pip install --editable .
```


## Usage
Run the autoredux command in the folder where the raw data lives
```sh
cd /path/to/raw/data
autoredux
```
