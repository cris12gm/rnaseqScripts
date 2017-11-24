#/usr/bin/env python3

import os,sys
import argparse
import math

##ARGPARSER

parser = argparse.ArgumentParser(description = 'Arguments')
parser.add_argument('--fcMatrix','-fc',type=str,help='FoldChange Matrix',required=True)
parser.add_argument('--filter','-filter',type=str,help='Comparations required with a |FC|>2',required=True)
parser.add_argument('--fpkm','-fpkm',type=str,help='FPKM matrix',required=True)
parser.add_argument('--go','-go',type=str,help='GO terms',required=True)
parser.add_argument('--outfile','-o',type=str,default='./output.filterFC',help='Output file')
parser.add_argument('--outfileG','-og',type=str,default='./outputGO.filterFC',help='Output file GOterms')
args = parser.parse_args()

##Abrir los ficheros y leer la cabecera

file_fc = open(args.fcMatrix)
file_filter = open(args.filter)
file_fpkm = open(args.fpkm)
file_go = open(args.go)

cabecera = file_fc.readline().strip().split("\t")
cabecera_fp = file_fpkm.readline().strip().split("\t")

fc = file_fc.readlines()
ff = file_filter.readlines()
fp = file_fpkm.readlines()
go = file_go.readlines()

##Creamos matrices y dict necesarios
toFilter = []
data = {}
fpkm = {}
goTerms = {}

##Guardamos en ff las comparaciones a filtrar(|FC|>2)
for line in ff:
	line = line.strip()
	toFilter.append(line)

##Guardamos los terminosGO
for line in go:
	line = line.strip().split("\t")
	for element in line:
		if "Ph" in element:
			gene = element
			goTerms[gene]=""
		else:
			goTerms[gene]=element

#Guardamos en fpkm
for line in fp:
	line = line.strip().split("\t")
	i = 0
	for element in line:
		if "Ph" in element:
			i = i+1
			gene = element
			m = {}
			fpkm[gene] = m
		else:
			sample = cabecera_fp[i]
			i = i+1
			m = fpkm[gene]
			m[sample] = element
			fpkm[gene] = m

##Guardamos en fc los datos de FoldChange
for line in fc:
	line = line.strip().split("\t")
	i = 0
	for element in line:
		if "Ph" in element:
			i = i+1
			gene = element
			m = {}
			data[gene]=m
		else:
			sample = cabecera[i]
			i = i+1
			m = data[gene]
			m[sample] = element
			data[gene] = m

results = [] #Matriz donde guardaré los genes que pasen mi filtro

##Recorremos la matriz de datos de FC

for gene in data:
	save = True #boolean que cambiaré a false en cuanto el gen no cumpla alguna condición
	m = data[gene] #en M están los datos de FC para ese gen
	for comparacion in m: #Miro todas los FC que hay en M
		if comparacion in toFilter: #Si esa comparación está en las del filtro
			value = abs(float(m[comparacion])) #guardo el valor de FC para esa comparacion en value
			if value<1:
				save = False #Si el gen en alguna de las comparaciones que tiene que tener |FC|>2 tiene menos, no se guarda
		else: #si la comparacion que estoy mirando no es de las |FC|<2
			value = abs(float(m[comparacion]))
			if value>1:
				save = False #Aquí queremos que |FC|<2, si no es, borro el gen
	if save == True: #si el gen ha pasado todos los filtros
		results.append(gene)

output = open(args.outfile,'a')
outputG = open(args.outfileG,'a')

for element in cabecera_fp:
	if "GENE" in element:
		linea_cabecera = str(element)
	else:
		linea_cabecera = linea_cabecera +"\t"+str(element)

linea_cabecera = linea_cabecera+"\n"
output.write(linea_cabecera)

for gene in results: #guardo los resultados en el fichero de salida
	linea = fpkm[gene]
	for element in cabecera_fp:
		if element=="GENE":
			linea_escribir = str(gene)
		else:
			linea_escribir = linea_escribir+"\t"+str(linea[element])
	linea_escribir = linea_escribir+"\n"
	output.write(linea_escribir)

	try:
		lineaGO = str(gene)+"\t"
		lineaGO = lineaGO+str(goTerms[gene])+"\n"
		outputG.write(lineaGO)
	except:
		pass
output.close()
outputG.close()
