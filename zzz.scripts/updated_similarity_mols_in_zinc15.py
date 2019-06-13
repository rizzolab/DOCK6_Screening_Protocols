### This script is passed a smiles string   ###
### searches ZINC using the url based       ###
### features, downloads the resulting html  ###
### and searches it for molecules.          ###


import sys
import urllib2
import urllib

### input
### command: python <this_script> <filename> <tan_cutoff> <zinc_return_count> <output filename>
filename = sys.argv[1]
tan_cutoff  = float(sys.argv[2]) #similarity threshold
zinc_return_count = int(sys.argv[3])
outputfilename = sys.argv[4]
print(filename,tan_cutoff,zinc_return_count)



file=open(filename, "r") 
lines=file.readlines()
file.close()

SMILES = []

### performs a logic check for spaces vs. commas a delimiters, based on the first line
### being 2 seperate columns
if len(lines[0].split(' ')) == 2:
   for line in lines:
      linesplit=line.split(' ')
      SMILES.append(linesplit[0].strip())
elif len(lines[0].split(',')) == 2:
   for line in lines:
      linesplit=line.split(',')
      SMILES.append(linesplit[0].strip())

else:
   print("Delimiters must be spaces or commas.")

file=open(outputfilename + ".dat" , "w")
for i in range(1,len(SMILES)):
   print (i)
   
   ### collected info
   matches = 0
   best_tan = 0
	

### http://zinc15.docking.org/substances/?highlight=CC%28C%29Cc1ccccc1&ecfp4_fp-tanimoto=CC%28C%29Cc1ccccc1
   pagename = "http://zinc15.docking.org/substances/?count=%d&ecfp4_fp-tanimoto-%f=%s" % (zinc_return_count,tan_cutoff,urllib.quote_plus(SMILES[i])) 
   print(pagename)
   try:
       page = urllib2.urlopen( pagename )
   except urllib2.HTTPError, error:
       print "empty", "empty", pagename
       continue
   j=0
   while(j==0): 
      for line in page:
           if '<div class="zinc-tile annotation small">' in line:
                   matches += 1 
           ### if (best_tan == 0):  #zinc rank orders starting with best tan, only need to save best one
           if "/substances/ZINC" in line:
                   line3 = line.replace("/" , ".")
                   spl3 = line3.split(".")
                   if (spl3[2] == "/n"):
                       continue
                   zinc_id = spl3[2]
                   print zinc_id
           if "<nobr>" in line:
                ### print line
                ### replace characters so we only have one character to split on
                   line2 = line.replace("<" , ">")
                   spl = line2.split(">")
                   if (spl[4] == "\n"):
                       continue
                   best_tan =  spl[4]
                   print best_tan
                   if float(best_tan)>=tan_cutoff:
                       print matches, best_tan, zinc_id, i, pagename
                       file.write('%s,%s,%s \n'% (SMILES[i],zinc_id,best_tan))
                   if float(best_tan)<tan_cutoff:
                       j=1
                   print j
           if j==1:
   	       break
file.close()
