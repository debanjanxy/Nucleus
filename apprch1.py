#!/usr/bin/python
#Nucleus
#Author: Debanjan Das		Author: Mamta Kumari
#Roll: 001310501011		Roll: 001310501009
#BCSE-II , A1			BCSE-II , A1

#JADAVPUR UNIVERSITY

#------------------------------------------
#START

import urllib2 
import urllib 
from bs4 import BeautifulSoup
import re
import string
import json
import time
import xlrd
from collections import OrderedDict
import operator



def getData(url):
    response = urllib.urlopen(url);
    time.sleep(1)
    data = json.loads(response.read())
    return data

#Searches using ensemble id in mamsap.it.deakin.edu.au and returns the responsible miRNAs as a list.
def Targets(ensembl):
	url = "http://mamsap.it.deakin.edu.au/~amitkuma/mirna_targetsnew/mirna_target_genes.php?gene_id_list=%0D%0A+++"+ensembl+"+&gene_name=ensembl_id&species=human&energy_cutoff=-30"
	data = urllib.urlopen(url)
	html = data.read()      
	txt = BeautifulSoup(html)
	data = txt.get_text()
	line = re.sub('[RM=]' , '' , data)
	List = line.split()
	List = [x for x in List if x != '--']
	word = List[2]
	word = word[23:35]
	del List[0]
	del List[1]
	List[0] = word
	return List

#Searches using gene name in genenames.org , gets the HGNC id, uses it to fetch the Ensemble ID.
#It goes to Targets() function and gets the list of miRNAs and returns the list to the main program.
def mi_RNA1(gene):
	url = "http://www.genenames.org/cgi-bin/search?search_type=symbols&search="+gene+"&submit=Submit"
	data = urllib.urlopen(url)
	html = data.read()      
	txt = BeautifulSoup(html)
	nice_html = txt.prettify()
	data = txt.get_text()
	line1 = re.sub('[.HGNC:]' , ' ' , data)
	list1 = line1.split()
	n = list1.index('_ID')
	id1 = list1[n+1]
	url = "http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:"+id1
	data = urllib.urlopen(url)
	html = data.read()
	txt = BeautifulSoup(html)
	data = txt.get_text()
	list2 = data.split()
	n1 = list2.index('Ensembl:')
	ensembl = list2[n1+1]
	miRNA1 = Targets(ensembl)
	return miRNA1



#Searches in TargetScan.org and returns the responsible miRNAs as a list
def mi_RNA2(gene):
	url = "http://www.targetscan.org/cgi-bin/targetscan/vert_61/targetscan.cgi?species=Human&gid=" + gene + "&mir_sc=&mir_c=&mir_nc=&mirg="
	data = urllib.urlopen(url)
	html = data.read()      
	txt = BeautifulSoup(html)
	data = txt.get_text()
	line1 = re.sub('[.,=!_@#${}:()&]' , ' ' , data)
	list1 = line1.split()
	n = list1.index('NM')
	num = list1[n + 1]
	str1 = "NM_" + num
	url = "http://www.targetscan.org/cgi-bin/targetscan/vert_61/view_gene.cgi?taxid=9606&rs=" + str1 + "&members=&showcnc=0&shownc=0&showncf="
	data1 = urllib.urlopen(url)
	html1 = data1.read()      
	txt1 = BeautifulSoup(html1)
	data1 = txt1.get_text()
	line2 = re.sub('[|.,=!_@#${}:()&]' , ' ' , data1)
	list2 = line2.split()
	indices = [i for i, x1 in enumerate(list2) if x1 == gene]
	for i1 in range(len(indices)):
		indices[i1] = indices[i1] + 2
	miRNA = list()
	for j in range(len(indices)):
		count = indices[j]
		miRNA.insert(j , list2[count])
	string0 = ' '.join(miRNA)
	string1 = re.sub('[UTR]' , '' , string0)
	miRNA = string1.split()
	return miRNA	

def Ensemble_ID(gene):
	url = "http://www.genenames.org/cgi-bin/search?search_type=symbols&search="+gene+"&submit=Submit"
	data = urllib.urlopen(url)
	html = data.read()      
	txt = BeautifulSoup(html)
	nice_html = txt.prettify()
	data = txt.get_text()
	line1 = re.sub('[.HGNC:]' , ' ' , data)
	list1 = line1.split()
	n = list1.index('_ID')
	id1 = list1[n+1]
	url = "http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:"+id1
	data = urllib.urlopen(url)
	html = data.read()
	txt = BeautifulSoup(html)
	data = txt.get_text()
	list2 = data.split()
	n1 = list2.index('Ensembl:')
	ensembl = list2[n1+1]
	return ensembl


def Atlas_Cancer(ensemble_id , gene_name):
	url = "http://www.proteinatlas.org/%s-%s/cancer" % (ensemble_id , gene_name)
	data = urllib.urlopen(url)
	html = data.read()      
	txt = BeautifulSoup(html)
	nice_html = txt.prettify()
	nice_html = re.sub('[ ]' , '' , nice_html)
	list1 = nice_html.split('>')
	n = len(list1)
	for i in range(len(list1)):
		if(list1[i] == '\n<tdclass="nowraptight"'):
			break
	list2 = list1[i : n+1]
	string1 = ''.join(list2)
	string1 = re.sub('[\<;>/]' , ' ' , string1)
	with open("atlas.txt" , "a") as myfile:
		myfile.write(string1)
	myfile = open('atlas.txt', 'r')
	text = myfile.read()
	myfile.close()
	f = open('atlas.txt', 'r+')
	f.truncate()
	text = re.sub('[\ \']+', " ", text)
	words = list(text.split())
	elem1 = 'tdclass="nowraptight"'
	elem2 = 'td'
	cancer_list = list()
	for i in range(len(words)):
		if(words[i] == elem1 and words[i+2] == elem2):
			cancer_list.append(words[i+1])
	cancer_list = list(OrderedDict.fromkeys(cancer_list))
	elem1 = '"onmouseover="tooltip('
	elem2 = ',0)'
	garbage_list = list()
	for i in range(len(words)):
		if(words[i] == elem1 and words[i+2] == elem2):
			garbage_list.append(words[i+1])
	garbage_string = '*'.join(garbage_list)
	str1 = garbage_string.replace(':' , '*')
	use_list = str1.split('*')
	score_list = list()
	for i in range(len(use_list)):
		if(use_list[i] == 'Staining'):
			score_list.append(use_list[i+1])
	n = len(cancer_list)
	score_list = score_list[0 : 4*n]
	str1 = ''.join(score_list)
	score_list = re.findall(r'\d+', str1)
	final_score = list()
	for i in range(len(score_list)):
		if( i%2 == 0 ):
			score = float(score_list[i]) / float(score_list[i+1])
			final_score.append(score)
	n = len(final_score)
	List = list()
	for i in range(0 , n , 4):
		first = 4 * final_score[i]
		second = 3 * final_score[i+1]
		third = 2 * final_score[i+2]
		fourth = 1 * final_score[i+3]
		sum = first + second + third + fourth
		List.append(sum)
	dict1 = dict(zip(cancer_list , List))
	sorted_x = sorted(dict1.items(), key=operator.itemgetter(1))
	n = len(sorted_x)
	sorted_x = sorted_x[n-5 : n]
	sorted_x = dict(sorted_x)
	final_cancer = list()
	for key in sorted_x.keys():
	   final_cancer.append(key)
	for i in range(len(final_cancer)):
		name = final_cancer[i]
		if(name != 'Lymphoma' and name != 'Carcinoid' and name != 'Melanoma' and name != 'Glioma'):
			n = len(name)
			name = name[0 : n-6]
			name = name + ' Cancer'
			final_cancer[i] = name
	return final_cancer

#Main Program
id = raw_input("Uniprot ID : ")

url = "http://mobidb.bio.unipd.it/ws/entries/%s/mobidb"%(id)
data = getData(url)
if len(data) == 1:
        print "ID doesnt exist."
else:
	print "ID exists as a Disordered Protein."
	data = urllib2.urlopen('http://www.uniprot.org/uniprot/'+ id +'.txt').read(1000)
	list1 = re.findall(r"[\w']+" , data)
	n = list1.index('Name')
	gene = list1[n + 1]
	print "\nGene Name : " + gene + "\n"
	ensembl = Ensemble_ID(gene)
	print "Ensemble ID : " + ensembl	
	Target_scan = mi_RNA2(gene)
	miRNA_target = mi_RNA1(gene)
	final = list()
	for i in range(len(Target_scan)):
		elem = Target_scan[i]
		for j in range(len(miRNA_target)):
			if(Target_scan[i] == miRNA_target[j]):
				final.append(elem)
		
	str1 = '*'.join(final)
	s = list(str1)
	ch = 'i'
	indices = [i for i, x1 in enumerate(s) if x1 == ch]
	for i in range(len(indices)):
		k = i + 1
		x = indices[i] + k
		s.insert(x , 'r')
	str1 = ''.join(s)
	final = str1.split('*')
	#for i in range(len(final)):
	#	print final[i]
	final = list(OrderedDict.fromkeys(final))
	cancer = list()
	atlas_cancer = Atlas_Cancer(ensembl , gene)
	wb = xlrd.open_workbook('alldata.xlsx')
	sh = wb.sheet_by_index(0)
	first_col = sh.col_values(0)
	second_col = sh.col_values(1)
	third_col = sh.col_values(2)
	cancer_result = list()
	mirna_result = list()
	for i in range(len(final)): 
		inp = final[i]
		indices1 = [i for i, x2 in enumerate(second_col) if x2 == inp]
		for j in range(len(indices1)):
			cnt = indices1[j]
			elem = third_col[cnt]
			for k in range(len(atlas_cancer)):
				if(atlas_cancer[k] == elem):
					mirna_result.append(inp)
					cancer_result.append(elem)
			cancer.append(elem)
	#cancer = list(OrderedDict.fromkeys(cancer))
	if(len(mirna_result) == 0 or len(cancer_result) == 0):
		print gene + " is not responsible for any kind of cancer.\n"
	else:
		final_dict = dict(zip(mirna_result , cancer_result))
		print "List of miRNAs vs Cancers : \n"
		for i in final_dict:
			print i + " : " + final_dict[i] + "\n"
	#if(len(cancer) != 0):
	#	print "\nProbable Cancers are :\n"
	#	for i in range(len(cancer)):
	#		print "\t" , cancer[i]
	#if(len(cancer) == 0):
	#	print "\nThere are no cancer diesease for this protein.\n"

		
	
	
	
