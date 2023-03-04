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

El archivo de configuracion tiene que estar en la carpeta configuration_files/ y debe ser de extension ".ini".

Ahi, se pueden especificar valores default y valores por ejercicios. Los valores de cada ejercicio deben ir debajo de su header adecuado ([ex1], [ex2], etc.).

No olvidarse de agregar todos los parametros, ya sea en el default o debajo del header del ejercicio. 

Para mas informacion de que parametros requiere cada ejercicio, puede ver el informe o correr el siguiente comando de ayuda:
```
python main.py --help
```
Los keys de los parametros en el archivo de configuracion deben tener el mismo nombre que los parametros al correrlo en linea de comandos (ej: --PARAM entonces el key es PARAM).

El comando a correr una vez que este todo seteado:
```
python main.py --exercise <num> --use_config --config_file <filename>
``` 
Donde ``` num``` es un numero de 1 a 5, ```use_config``` indica que se quiere usar un archivo de configuracion y ```config_file``` es que archivo se quiere usar, el cual es ```configuration.ini``` por default.
