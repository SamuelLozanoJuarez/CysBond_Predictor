
#importamos los paquetes necesarios
#si alguno no está instalado lo instala
import os
import argparse
import warnings
warnings.filterwarnings("ignore")

#importamos Biopython
try:
    import Bio
except ImportError:
    print("The Biopython package is required. It will be installed...")
    try:
        import pip
        pip.main(['install', 'biopython'])
        import Bio
        print("Biopython package successfully installed")
    except Exception as e:
        print("An error occurred during the installation of Biopython:", e)

#importamos Numpy
try:
    import numpy as np
except ImportError:
    print("The Numpy package is required. It will be installed...")
    try:
        import pip
        pip.main(['install', 'numpy'])
        import numpy as np
        print("Numpy package successfully installed")
    except Exception as e:
        print("An error occurred during the installation of Numpy:", e)

#importamos statistics
try:
    import statistics as stats
except ImportError:
    print("The statistics package is required. It will be installed...")
    try:
        import pip
        pip.main(['install', 'statistics'])
        import statistics as stats
        print("statistics package successfully installed")
    except Exception as e:
        print("An error occurred during the installation of statistics:", e)

#importamos tabulate
try:
    import tabulate
except ImportError:
    print("The tabulate package is required. It will be installed...")
    try:
        import pip
        pip.main(['install', 'tabulate'])
        import tabulate
        print("tabulate package successfully installed")
    except Exception as e:
        print("An error occurred during the installation of tabulate:", e)

#cargamos las funciones necesarias
from tabulate import tabulate
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral


#leemos los parámetros dados en la ejecución del script
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True, help='Relative or absolute path to the PDB file')
parser.add_argument('-o', '--output', type=str, default='', help='Relative path where the PyMOL script is stored (by default: current directory)')
parser.add_argument('-a', '--angle', type=str, default='(84,96)', help='Minimum and maximum value of the dihedral angle between cysteines to consider disulphide bridging potential. Must be enclosed in parentheses and separated by a comma.')
parser.add_argument('-d', '--distance', type=str, default='(1.5,2.5)', help='Minimum and maximum value of the distance between cysteines to consider disulphide bridging potential. Must be enclosed in parentheses and separated by a comma.')
parser.add_argument('-aci', '--angle_ci', type=int, default=5, help='Value of the confidence interval to be considered above and below the maximum and minimum values of the dihedral angle. Must be an integer.')
parser.add_argument('-dci', '--distance_ci', type=float, default=0.1, help='Value of the confidence interval to be considered above and below the maximum and minimum values of the distance. Must be a float.')
parser.add_argument('-sci', '--show_ci', type=str, default='False', help='Boolean variable. Indicates whether to display the confidence interval bonds in PyMOL or not.')
#parseamos los argumentos
args = parser.parse_args()

#cargamos el valor de los argumentos en el script
file = args.input
output_path = args.output
angle_inf = int(args.angle.split(',')[0].replace('(',''))
angle_sup = int(args.angle.split(',')[1].replace(')',''))
distance_inf = float(args.distance.split(',')[0].replace('(',''))
distance_sup = float(args.distance.split(',')[1].replace(')',''))
low_ci_angle = angle_inf - args.angle_ci
high_ci_angle = angle_sup + args.angle_ci
low_ci_distance = distance_inf - args.distance_ci
high_ci_distance = distance_sup + args.distance_ci
show_ci = True if args.show_ci == 'True' else False


#construimos la ruta absoluta hasta el archivo si es necesario
if ':' not in file:
    abs_path = os.getcwd().replace('\\','/') + '/' + file.replace('\\','/')
else:
    abs_path = file

#cargamos el archivo PDB 
parser = PDBParser()
structure = parser.get_structure("",abs_path)
#comprobamos si es una estructura obtenida experimentalmente o predicha, a partir de la información del encabezado del PDB 
if 'alphafold' in structure.header['name'].lower():
    is_experimental = False
else:
    is_experimental = True

#definimos una función que dada una estructura de PDB devuelva una lista con todos los resiudos de cisteína
def get_cys(structure):
    '''
    Dado un objeto de tipo Structure devuelve todos los residuos de Cisteína que encuentre en él.

    Parameters
    ---------
     - structure: instancia de la clase Bio.PDB.Structure.structure que contiene la estructura descrita en el archivo PDB que queremos analizar.

    Return
    ---------
     Devuelve una lista de objetos de tipo Bio.PDB.Residue.residue que constituyen todos los residuos de cisteína encontrados en la estructura pasada como parámetro.
    '''
    #inicializamos una lista vacía para almacenar las cisteínas
    cys_list = []
    #recorremos las cadenas del archivo PDB
    for chain in structure.get_chains():
        #recorremos los residuos de la cadena
        for residue in chain.get_residues():
            #comprobamos si su nombre coincide con el de una cisteína
            if residue.get_resname() == 'CYS':
                #si es así los añadimos a la lista de cisteínas
                cys_list.append(residue)
    #devolvemos la lista
    return cys_list

#ahora definimos una función que dada una lista de residuos de cisteína y si es experimental o predicho, y elimina aquellas potencialmente móviles
def filter_cys(cys_list, is_experimental):
    '''
    Dada una lista con residuos de cisteína elimina aquellos potencialmente móviles según su B-factor o pLDDT

    Parameters
    ---------
     - cys_list: lista de residuos de cisteína (instancias de Bio.PDB.Residue.residue) que se quiere filtrar
     - is_experimental: booleano que indica si la estructura es experimental o no (predicha)
    '''
    if is_experimental:
        #si la estructura es experimental eliminamos aquellos residuos cuya media de B-factor sea superior a 30
        cys_filtered = list(filter(lambda cys: stats.mean(map(lambda atom: atom.bfactor, cys.get_atoms()))<=30, cys_list))
    else:
        #si la estructura no es experimental eliminamos los residuos con valor pLDDT inferior a 50
        cys_filtered = list(filter(lambda cys: stats.mean(map(lambda atom: atom.bfactor, cys.get_atoms()))>=50, cys_list))
    #devolvemos la lista filtrada
    return cys_filtered

#definimos una función que devuelve el ángulo diedro entre los dos planos formados por los átomos de las cisteínas
def dihedral_cys(cys1, cys2):
    '''
    Dados dos residuos de cisteína devuelve el ángulo diedro formado entre los átomos de C-Beta y S de ambos residuos.

    Parameters
    ---------
     - cys1: instancia de Bio.PDB.Residue.residue que representa el primer residuo de cisteína considerado para el ángulo.
     - cys2: instancia de Bio.PDB.Residue.residue que representa el segundo residuo de cisteína considerado para el ángulo.

    Return
    ---------
     Devuelve un número decimal positivo que se corresponde con los grados del ángulo diedro entre los dos planos formados por los residuos.
    '''
    #primero obtenemos los átomos de interés de cada cisteína (CB se corresponde con Carbono Beta y SG con el átomo de azufre)
    CB1, SG1 = list(filter(lambda atom: atom.id == 'CB' or atom.id == 'SG',cys1.get_unpacked_list()))
    CB2, SG2 = list(filter(lambda atom: atom.id == 'CB' or atom.id == 'SG',cys2.get_unpacked_list()))
    #y empleamos la función calc_dihedral que incorpora PDB para calcular el ángulo (en valor absoluto y en grados)
    return abs(np.degrees(calc_dihedral(CB1.get_vector(), SG1.get_vector(), SG2.get_vector(), CB2.get_vector())))

#haciendo uso de las funciones anteriores obtenemos una lista de cisteínas que cumplen los requisitos de calidad (B-factor y pLDDT)
filtered_cys = filter_cys(get_cys(structure), is_experimental)

#iniciamos una lista para añadir los pares de cisteínas que cumplan los requisitos de distancia y ángulo para formar potencial puente diS 
bonds = []
#y una lista para almacenar los pares de cisteínas que no cumplan los requisitos estrictos pero sí el intervalo de confianza
ci_bonds = []
#recorremos la lista de cisteínas ya filtradas
for cys1 in filtered_cys:
    #obtenemos el átomo de S correspondiente a la cisteína de esta iteración
    SG_atom1 = list(filter(lambda atom: atom.get_name() == 'SG', cys1.get_unpacked_list()))[0]
    #recorremos internamente otra vez la lista de cisteínas
    for cys2 in filtered_cys:
        #obtenemos el átomo de S de la cisteína de la segunda iteración
        SG_atom2 = list(filter(lambda atom: atom.get_name() == 'SG', cys2.get_unpacked_list()))[0]
        #comprobamos si la distancia y el ángulo diedro están en el rango de valores deseados
        if distance_inf<=(SG_atom1 - SG_atom2)<=distance_sup and angle_inf<=dihedral_cys(cys1,cys2)<=angle_sup:
            #si es así vamos a crear un diccionario con la información del potencial puente diS
            info = dict()
            #añadimos información de la primera cisteína
            info['Cys1'] = {'Name':cys1.get_resname(), 'Pos':cys1.id[1]}
            #de la segunda
            info['Cys2'] = {'Name':cys2.get_resname(), 'Pos':cys2.id[1]}
            #y el valor de distancia y ángulo
            info['Characteristics'] = {'Distance':(SG_atom1 - SG_atom2), 'Dihedral Angle':dihedral_cys(cys1,cys2)}
            bonds.append(info)
        elif low_ci_distance<=(SG_atom1 - SG_atom2)<=high_ci_distance and low_ci_angle<=dihedral_cys(cys1,cys2)<=high_ci_angle:
            #si es así vamos a crear un diccionario con la información del potencial puente diS
            info = dict()
            #añadimos información de la primera cisteína
            info['Cys1'] = {'Name':cys1.get_resname(), 'Pos':cys1.id[1]}
            #de la segunda
            info['Cys2'] = {'Name':cys2.get_resname(), 'Pos':cys2.id[1]}
            #y el valor de distancia y ángulo
            info['Characteristics'] = {'Distance':(SG_atom1 - SG_atom2), 'Dihedral Angle':dihedral_cys(cys1,cys2)}
            ci_bonds.append(info)


#a continuación debemos eliminar los duplicados de la lista, ya que un mismo enlace CYS1-CYS2 se ha considerado también como CYS2-CYS1
unique_combinations = set()
#inicializamos una nueva lista para almacenar los puentes únicos
unique_bonds = []
#recorremos la lista de enlaces diS
for item in bonds:
    #obtenemos la pareja de cisteínas entre las que se forma potencialmente el enlace
    cys1_pos = item['Cys1']['Pos']
    cys2_pos = item['Cys2']['Pos']
    combination = tuple(sorted((cys1_pos, cys2_pos)))
    #comprobamos si este par de cisteínas ya está contemplado
    if combination not in unique_combinations:
        #si no ha sido contemplado hasta ahora actualizamos la lista de puentes y la de pares de cisteínas contempladas
        unique_combinations.add(combination)
        unique_bonds.append(item)

#hacemos el mismo proceso con los posibles enlaces que hayan caído en el intervalo de confianza
unique_combinations_ci = set()
#inicializamos una nueva lista para almacenar los puentes únicos
unique_bonds_ci = []
#recorremos la lista de enlaces diS
for item in ci_bonds:
    #obtenemos la pareja de cisteínas entre las que se forma potencialmente el enlace
    cys1_pos = item['Cys1']['Pos']
    cys2_pos = item['Cys2']['Pos']
    combination = tuple(sorted((cys1_pos, cys2_pos)))
    #comprobamos si este par de cisteínas ya está contemplado
    if combination not in unique_combinations_ci:
        #si no ha sido contemplado hasta ahora actualizamos la lista de puentes y la de pares de cisteínas contempladas
        unique_combinations_ci.add(combination)
        unique_bonds_ci.append(item)

#y ya solo queda mostrar la información por pantalla
#vamos a usar el paquete tabulate para hacerlo todo más elegante
headers = ['Bond','Cys1','Cys2','Distance','Dihedral Angle']
out_table = []
for i, item in enumerate(unique_bonds):
    #creamos una lista con los valores de la fila correspondiente a este puente diS (cisteínas, distancia, ángulo)
    cys1 = ' '.join([str(i) for i in item['Cys1'].values()])
    cys2 = ' '.join([str(i) for i in item['Cys2'].values()])
    dist = item['Characteristics']['Distance']
    angle = item['Characteristics']['Dihedral Angle']
    #y lo añadimos a la tabla que vamos a mostrar
    out_table.append([i+1 ,cys1, cys2, dist, angle])

#también creamos la tabla para los enlaces del intervalo de confianza
out_table_ci = []
for i, item in enumerate(unique_bonds_ci):
    #creamos una lista con los valores de la fila correspondiente a este puente diS (cisteínas, distancia, ángulo)
    cys1 = ' '.join([str(i) for i in item['Cys1'].values()])
    cys2 = ' '.join([str(i) for i in item['Cys2'].values()])
    dist = item['Characteristics']['Distance']
    angle = item['Characteristics']['Dihedral Angle']
    #y lo añadimos a la tabla que vamos a mostrar
    out_table_ci.append([i+1 ,cys1, cys2, dist, angle])

#finalmente comprobamos si hay algún potencial puente diS que mostrar
if out_table:
    print('\nThe following potential disulfide bonds have been detected:')
    print(tabulate(out_table, headers=headers, tablefmt = "pretty"))
#si no hay ningún puente diS se muestra un mensaje por pantalla advirtiendo al usuario
else:
    print('\n --------------------------------- ')
    print('No potential disulfide bond has been detected.')
    print(' --------------------------------- ')

#y lo mismo con los enlaces que estén en el intervalo de confianza
if out_table_ci:
    print('\nThe following disulfide bonds do not meet the requirements but should be reviewed:')
    print(tabulate(out_table_ci, headers=headers, tablefmt = "pretty"))
#si no hay ningún puente diS se muestra un mensaje por pantalla advirtiendo al usuario
else:
    print('\n --------------------------------- ')
    print('There are no additional disulphide bonds that need to be checked.')
    print(' --------------------------------- ')

#a modo de añadido vamos a generar un archivo .py que pueda ser ejecutado en PyMOL para visualizar los puentes disulfuro encontrados, así como sus características
#inicializamos una lista para almacenar las líneas del script
script_lines = [
    '#primero importamos pymol',
    'from pymol import cmd\n',
    '#cargamos el archivo PDB analizado',
    f'cmd.load("{abs_path}")'
]


#antes de generar los scripts debemos comprobar la variable show_ci, para determinar si deben incluirse en PyMOL los enlaces que hayan caído en el intervalo de confianza
if show_ci:
    unique_bonds += unique_bonds_ci


#recorremos los puentes disulfuro
for bond in unique_bonds:
    #para cada puente obtenemos la posición de ambas cisteínas
    pos1 = bond['Cys1']['Pos']
    pos2 = bond['Cys2']['Pos']
    script_lines.append(f'\n#Puente disulfuro {list(bond["Cys1"].values())} - {list(bond["Cys2"].values())}')
    script_lines.append('\n#seleccionamos los atomos de interes de cada cisteina')
    #y añadimos a la lista el comando para seleccionar cada uno de los 3 átomos de interés de cada cisteína (C-Beta, C-Alfa y S)
    for position in [pos1,pos2]:
        script_lines.append(f"cmd.select('SG{position}', 'resi {position} and name SG')")
        script_lines.append(f"cmd.select('CB{position}', 'resi {position} and name CB')")
        script_lines.append(f"cmd.select('CA{position}', 'resi {position} and name CA')")
    script_lines.append('\n#creamos los enlaces entre atomos')
    #y añadimos también el comando para crear los enlaces cBeta-cAlfa, cBeta-azufre y azufre-azufre y mostrarlos
    script_lines.append(f"cmd.bond('CA{pos1}','CB{pos1}')")
    script_lines.append(f"cmd.show('sticks','CA{pos1} or CB{pos1}')")
    script_lines.append(f"cmd.bond('CB{pos1}','SG{pos1}')")
    script_lines.append(f"cmd.show('sticks','CB{pos1} or SG{pos1}')")
    script_lines.append(f"cmd.bond('SG{pos1}','SG{pos2}')")
    script_lines.append(f"cmd.show('sticks','SG{pos1} or SG{pos2}')")
    script_lines.append(f"cmd.bond('SG{pos2}','CB{pos2}')")
    script_lines.append(f"cmd.show('sticks','SG{pos2} or CB{pos2}')")
    script_lines.append(f"cmd.bond('CB{pos2}','CA{pos2}')")
    script_lines.append(f"cmd.show('sticks','CB{pos2} or CA{pos2}')")
    #añadimos el comando para mostrar la distancia entre las cisteínas
    script_lines.append("\n#calculamos la distancia entre cisteinas")
    script_lines.append(f"cmd.distance('','SG{pos1}','SG{pos2}')")
    #y el comando para mostrar el ángulo diedro entre los dos planos
    script_lines.append("#y el angulo diedro entre los dos planos")
    script_lines.append(f"cmd.dihedral('','CB{pos1}','SG{pos1}','SG{pos2}','CB{pos2}')")
    #y el nombre de los residuos
    script_lines.append("#y el nombre de las cisteinas involucradas")
    script_lines.append(f"cmd.label('n. SG and i. {pos1}', '\"SG-CYS%s\" %resi')")
    script_lines.append(f"cmd.label('n. SG and i. {pos2}', '\"SG-CYS%s\" %resi')")


#posteriormente, si se ha encontrado algún puente disulfuro se añade un comando para hacer zoom en el primer enlace encontrado
if unique_bonds:
    script_lines.append("#finalmente hacemos zoom en el primer enlace")
    pos = unique_bonds[0]['Cys1']['Pos']
    script_lines.append(f"cmd.zoom('SG{pos}')")


#comprobamos si la ruta de salida tiene el aspecto requerido, sino lo modificamos para que sea así
#básicamente comprueba que si se ha proporcionado ruta de salida esta sea un directorio terminado en '/'
if output_path!='' and output_path[-1] != '/':
    output_path += '/'
    #y si no existe el directorio de salida lo crea
    if not os.path.exists(output_path):
        os.mkdir(output_path)

#y finalmente escribimos el script en el mismo directorio desde el que se ejecuta el código o en su defecto en el parámetro -o 
with open(f'{output_path}{abs_path.split("/")[-1].split(".")[0]}.py','w') as archivo:
    for linea in script_lines:
        archivo.write(f"{linea}\n")

