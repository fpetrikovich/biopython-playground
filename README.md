# Biopython Playground

## Installation

Create a virtual environment with `virtualenv`:
```
virtualenv .env
```
or
```
python3 -m venv "virtualenv"
```

Activate the environment:
```
source .env/bin/activate
```

Install required dependencies:
```
pip install -r requirements.txt
```

## Using the configuration file

The configurations files must be in the directory `configuration_files` and must have the .ini extension.

You can run the ORF analysis with the default configuration file (`configuration.ini`) as follows:
```
python main.py --analysis orf
```