# Padé-Rayligh-Ritz

## Auteurs

* **PAILLEUX Jean-Didier**
* **JOULIN Maxence**
* date: 27/02/2019 <br>

## Compilation
La compilation du projet se fait de la manière suivante:
```
make
```

## Execution
L'exécution du code PRR en parallèle en utilisant K processus MPI avec P itération sur une matrice N*N se fait de la manière suivante:
* Pour des matrice provenant de Matrix Market:
```
mpirun -n 2 ./$(EXEC) mm matrix1.mtx  
```
* Pour des matrices au format txt:
```
mpirun -n 2 ./$(EXEC) txt matrix2.txt  
```
