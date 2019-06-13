### This script is passed a smiles string   ###
### searches PUBCHEM using the url based    ###
### features, downloads the resulting html  ###
### and searches it for molecules.          ###


import sys
import urllib2
import urllib
import time
### input
### command: python <this_script> <filename> <tan_cutoff> <number_keep> <outputfilename>
filename = sys.argv[1] #input filename
tan_cutoff  = float(sys.argv[2]) #similarity threshold
number_keep = int(sys.argv[3]) #number of cid's to keep from pubchem
outputfilename = sys.argv[4] #output filename 


file=open(filename, "r") 
lines=file.readlines()
file.close()

SMILES = []
identifiers = []
### performs a logic check for spaces vs. commas a delimiters, based on the first line
### being 2 seperate columns
if len(lines[0].split(' ')) == 2:
   for line in lines:
      linesplit=line.split(' ')
      SMILES.append(linesplit[0].strip())
      identifiers.append(linesplit[1].strip())
elif len(lines[0].split(',')) == 2:
   for line in lines:
      linesplit=line.split(',')
      SMILES.append(linesplit[0].strip())
      identifiers.append(linesplit[1].strip())

else:
   print("Delimiters must be spaces or commas.")

file=open(outputfilename+".dat", "w")
for i in range(1,len(SMILES)):
   print (i)
   matches = 0

### http://zinc15.docking.org/substances/?highlight=CC%28C%29Cc1ccccc1&ecfp4_fp-tanimoto=CC%28C%29Cc1ccccc1
   pagename = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/%s/JSON?Threshold=%d" % (urllib.quote_plus(SMILES[i]), tan_cutoff)
   print(pagename)
   cid_scrape_list=[]
   smiles_scrape_list=[]
   try:
      page = urllib2.urlopen( pagename )
   except urllib2.HTTPError as error:
      error_message = error.read()
      if "PUGREST.Timeout" in error_message:
         file.write(identifiers[i]+","+SMILES[i]+","+"Timeout Error"+"\n")
      if "PUGREST.NotFound" in error_message:
         file.write(identifiers[i]+","+SMILES[i]+","+"No compounds found for cutoff."+"\n")
      continue
   ### This line adds the de novo compound line - can comment out if necessary
   file.write(identifiers[i]+","+SMILES[i]+","+"DE_NOVO_COMPOUND"+"\n")
   pagelist=[]
   for line in page:
      pagelist.append(line)
   for i in range(len(pagelist)):
      if "cid" in pagelist[i]:
         splitline=pagelist[i].split(":")
         cid_scrape_list.append(splitline[1].strip())
         matches+=1
      if ("SMILES" in pagelist[i]) and ("Canonical" in pagelist[i+1]):
         splitline=pagelist[i+9].split(":")
         smiles_scrape_list.append(splitline[1].replace('"','').strip())
      if len(smiles_scrape_list) == number_keep:
         break
   combined_list=[]
   for i in range(0,len(cid_scrape_list)):
      combined_list.append([cid_scrape_list[i],smiles_scrape_list[i]])
   for pair in combined_list:
      file.write(str(pair[0])+","+str(pair[1]+"\n"))
   time.sleep(0.20)
file.close()
