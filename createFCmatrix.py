#usr/bin/env python3

import os,sys
import argparse
import math

####ARGPARSER

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--fpkmMatrix','-f',type=str,help='FPKM matrix in xls format',required=True)
parser.add_argument('--comparations','-c',type=str,help='Comparations file',required=True)
parser.add_argument('--outfile','-o',type=str,default='./outputFC.expr',help='Output file')
args=parser.parse_args()

f_fpkm = open(args.fpkmMatrix)
f_comp = open(args.comparations)

cabecera = f_fpkm.readline().strip().split("\t")

f = f_fpkm.readlines()
c = f_comp.readlines()

compare = []
data = {}

for line in c:
	line = line.strip()
	compare.append(line)

for line in f:
	line = line.strip().split("\t")
	i = 0
	for element in line:
		if "Ph" in element:
			gene = element
			m = {}
			data[gene] = m
		else:
			sample = cabecera[i]
#			sample = cabecera[i].split(".bam")[0]
			i = i+1
			m = data[gene]
			m[sample] = element
			data[gene] = m

results = {}

for gene in data:
	m_result= {}
	results[gene]=m_result
	m = data[gene]
	for comparation in compare:
		sample_a = comparation.split(",")[0]
		sample_b = comparation.split(",")[1]

		fpkm_a = float(m[sample_a])+1
		fpkm_b = float(m[sample_b])+1
		cociente = float(fpkm_a/fpkm_b)

		log2_foldChange = math.log(cociente,2)
		m_result = results[gene]
		m_result[comparation]=log2_foldChange

	results[gene] = m_result

cabecera_output=[]
##Create output first line

output = open(args.outfile,"a")
linea_c = "GENE"

for key in m_result:
	cabecera_output.append(key)
	linea_c = linea_c+"\t"+key

linea_c = linea_c+"\n"
output.write(linea_c)

for key in results:
	linea_o=""
	m_results = results[key]
	linea_o = str(key)
	for c in cabecera_output:
		resultado = m_results[c]
		linea_o = linea_o+"\t"+str(resultado)

	linea_o = linea_o+"\n"
	output.write(linea_o)
