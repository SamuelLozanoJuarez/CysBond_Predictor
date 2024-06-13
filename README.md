![logo](https://github.com/SamuelLozanoJuarez/CysBond_Predictor/assets/80346399/6b7513ed-aba8-4028-a2c2-8db3498c9ba6)


# CysBond Predictor

## ¿Qué es CysBond Predictor?

Se trata de un script de Python que permite identificar potenciales puentes disulfuro en un archivo PDB que contenga la estructura de una proteína. Además de reportar estos puentes disulfuro también genera un fichero .py que permite la visualización de la estructura tridimensional de la proteína y los puentes disulfuro utilizando la herramienta PyMOL.
	El programa permite al usuario personalizar los parámetros que determinan los potenciales puentes disulfuro (distancia entre cisteínas y ángulo diedro formado). Para que un par de cisteínas sean consideradas para un potencial puente disulfuro deben cumplir por defecto con un ángulo de entre 84º y 96º y una distancia de entre 1.5 Å y 2.5 Å. Por motivos de calidad de la estructura, las cisteínas con un B-factor superior a 30 Å $^2$ si la estructura se ha obtenido experimentalmente, o un pLDDT inferior a 50 si la estructura ha obtenido por predicción, serán descartadas.

## ¿Qué se incluye en el repositorio?


## Dependencias necesarias

Para poder ejecutar el script es necesario tener instalados los siguientes programas y paquetes. Se incluyen las versiones utilizadas en el desarrollo del script, con las que no hay problemas de compatibilidades:
-	_Python_ versión 3.11.8
-	_PyMOL_ versión 3.0.2
-	Paquete _Biopython_ versión 1.83
-	Paquete _NumPy_ versión 1.26.4
-	Paquete _tabulate_ versión 0.9.0

## Contacto y mantenimiento

Este programa ha sido desarrollado por Samuel Lozano Juárez, como parte de la asignatura Bioinformática Estructural del Máster en Bioinformática de la Universitat de València (España). Para cualquier consulta acerca del funcionamiento del programa dirigirse a una de las siguientes direcciones de correo: 

salojua@alumni.uv.es | samuellozanojuarez@gmail.com

El software es de libre uso para quien lo desee.
