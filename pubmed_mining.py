import csv
genes=[];

with open('/Users/andrejeremcuk/Downloads/genes.txt', 'r') as fp :
    reader = csv.reader(fp, delimiter='\t')
    for i in range(20000): 
     genes.append(reader.next())

genesq=[];
for i in range(len(genes)):
 genesq.append(genes[i][1])

del genesq[0]
del genes[0]
np.savetxt('/Users/andrejeremcuk/Downloads/genesq.txt', genesq,fmt='%s')

#######################################################################################
#######################################################################################

import time
import numpy as np
genesq=np.genfromtxt('/Users/andrejeremcuk/Downloads/genesq.txt',dtype='str')

from Bio import Entrez
from Bio import Medline

MAX_COUNT = 100
Entrez.email = 'a-nai@yandex.ru'
articles=[];genes_cancer_poor=[];genes_cancer_poor1=[];


for u in range(0,len(genesq)):
 print u
 if u%100==0: 
  np.savetxt('/Users/andrejeremcuk/Downloads/genes_cancer_poor.txt', genes_cancer_poor,fmt='%s');
  np.savetxt('/Users/andrejeremcuk/Downloads/genes_cancer_poor1.txt', genes_cancer_poor1, fmt='%s')
 gene=genesq[u];genefullname=genes[u][2]
 TERM=gene+' '+'poor prognosis'
 try: h=Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
 except: time.sleep(5);h=Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
 result = Entrez.read(h)
 ids = result['IdList']
 h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
 ret = Medline.parse(h)
 fer=[];
 for re in ret:
  try: tr=re['TI'];
  except: tr='0';
  fer.append(tr);
 for i in range(len(fer)):
  gene1=fer[i].find(gene)
  gene2=fer[i].find(genefullname)
  #####
  inc=fer[i].find("Increased")
  highe=fer[i].find("High expression")
  high=fer[i].find("High")
  expr=fer[i].find("expression")
  Overe=fer[i].find("Overexpression")
  overe=fer[i].find("overexpression")
  up1=fer[i].find("Up-regulation")
  el1=fer[i].find("Elevated expression")
  expr1=fer[i].find("Expression of ")
  ####
  decr=fer[i].find("Decreased")
  loss=fer[i].find("Loss")
  low1=fer[i].find("Low expression")
  low2=fer[i].find("Low levels")
  down1=fer[i].find("Down-regulated")
  down2=fer[i].find("Down-regulated")
  down3=fer[i].find("Downregulation")
  #####
  acc=fer[i].find("accelerates")
  poor=fer[i].find("poor patient prognosis")
  poor1=fer[i].find("poor prognosis")
  poor2=fer[i].find("unfavorable clinical outcomes")
  poor3=fer[i].find("unfavorable prognosis")
  poor4=fer[i].find("poor outcome")
  poor5=fer[i].find("poor survival")
  poor6=fer[i].find("poor patient survival")
  poor7=fer[i].find("progression and prognosis")
  ###
  canc=fer[i].find("cancer")
  canc1=fer[i].find("carcinoma")
  ###############################################
  if (gene1!=-1)or(gene2!=-1): #<poor1,poor,poor2,poor3,poor4,poor5,poor6,poor7
   if (canc1!=-1)or(canc!=-1):
    if (poor!=-1)or(poor1!=-1)or(poor2!=-1)or(poor3!=-1)or(poor4!=-1)or(poor5!=-1)or(poor6!=-1)or(poor7!=-1): #
     genel=-1;
     if (gene1!=-1): genel=gene1;
     if (gene2!=-1): genel=gene2;
     gene1=genel;
     if (expr!=-1): #<poor1,poor,poor2,poor3,poor4,poor5,poor6,poor7
      if (gene1<expr): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,1))
     if (low1!=-1)and(gene1!=-1):
      if (low1<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,2))
     if (el1!=-1)and(gene1!=-1):
      if (el1<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,3))
     if (Overe!=-1)and(gene1!=-1):
      if (Overe<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,4))
     if (overe!=-1)and(gene1!=-1):
      articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,5))
     if (expr1!=-1)and(gene1!=-1):
      if (expr1<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,6))
     if (up1!=-1)and(gene1!=-1):
      if (up1<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,7))
     if (highe!=-1)and(gene1!=-1):
      if (highe<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,8))
     if (high!=-1)and(gene1!=-1)and(expr!=-1):
      if (high<gene1<expr): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,9))
     if (gene1!=-1)and(expr1!=-1):
      if (expr1<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,10))
     if (gene1!=-1)and(inc!=-1):
      if (inc<gene1): 
       articles.append((fer[i],TERM,gene,u,i));genes_cancer_poor.append((gene,u,i,11))
     ###########
     if (gene1!=-1)and(decr!=-1):
      if (decr<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,12))
     if (gene1!=-1)and(loss!=-1):
      if (loss<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,13))
     if (gene1!=-1)and(low1!=-1):
      if (low1<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,14))
     if (gene1!=-1)and(low2!=-1):
      if (low2<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,15))
     if (gene1!=-1)and(down1!=-1):
      if (down1<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,16))
     if (gene1!=-1)and(down2!=-1):
      if (down2<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,17))
     if (gene1!=-1)and(down3!=-1):
      if (down3<gene1): 
       articles.append((fer[i],TERM,gene,u,i,'low'));genes_cancer_poor1.append((gene,u,i,18))


np.savetxt('/Users/andrejeremcuk/Downloads/articles.txt', articles,fmt='%s', delimiter='<');

