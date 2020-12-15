#!/usr/bin/env python
# coding: utf-8

# # This script needs Mitty and gargammel to be installed and running properly, exported from Jupyter Notebook

# ### Import packages

# In[1]:


import subprocess
import os
import random


# ### Set all variables

# In[2]:


ReadsWanted=1500000 #First we simulate reads from 20 to 170 (15 bins) bp but only use a subset for our analysis
do_linear_simulation = True
coverage=4
simulation_vcf_file = "1kg.22.NA19471_Indel.vcf.gz" #from "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" and we kept only the wanted individuals, in this case the African individual NA19471
linear_ref_fasta="human_g1k_v37.fa.gz" #from ftp "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
SampleName="NA19471"
Prefix="NA19471-unpaired-reads"
PathGargammel="Software/gargammel"



# ### Cutting reads into the length bin we want to have the same numbers in each bins. 

# In[10]:


def cutting_modern_reads_to_aDNA(FileInp,FileOut,Nb_reads):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator       
    #15 bins, compute the amount of reads needed in each bins
	RLengthTrue = {
        "20-29" : Nb_reads/15.0,
        "30-39" : Nb_reads/15.0,
        "40-49" : Nb_reads/15.0,
        "50-59" : Nb_reads/15.0,
        "60-69" : Nb_reads/15.0,
        "70-79" : Nb_reads/15.0,
        "80-89" : Nb_reads/15.0,
        "90-99" : Nb_reads/15.0,
        "100-109" : Nb_reads/15.0,
        "110-119" : Nb_reads/15.0,
        "120-129" : Nb_reads/15.0,
        "130-139" : Nb_reads/15.0,
        "140-149" : Nb_reads/15.0,
        "150-159" : Nb_reads/15.0,
        "160-169" : Nb_reads/15.0
    }
    
    if(sum(RLengthTrue.values()) != Nb_reads):
        print("ERROR ! Issue with your % !")
    else:
        print("Empirical distribution created")
        
    #Change the distribution to a % corresponding to our input file
    LineCount =subprocess.check_output(['wc', '-l', FileInp])
    LineCount=list(filter(None,LineCount.decode("utf-8").split(" ")))
    LineCount=int(LineCount[0])/int(4)
    for key in RLengthTrue:
        RLengthTrue[key] = round(RLengthTrue[key] * LineCount / Nb_reads)

    if(sum(RLengthTrue.values()) != LineCount):
        print("ERROR ! Probably a round problem !")
    else:
        print("Distribution transformed correctly")
    
    #Picking randomly a bin and a value between the extreme value within the bin to cut the reads (and the Q values) in the corresponding length
    handle = open(FileOut, "w")
    for title, seq, qual in FastqGeneralIterator(open(FileInp)) :
        #randomly pick a length in the dict
        length, number = random.choice(list(RLengthTrue.items()))
        while(number==0):
            del RLengthTrue[length]
            length, number = random.choice(list(RLengthTrue.items()))
        #Updating the dict 
        RLengthTrue[length]=RLengthTrue[length]-1
        if(RLengthTrue[length]==0):
            del RLengthTrue[length]
        #Use the next line only if it is a range in the distribution
        length=random.randint(int(length.split("-")[0]), int(length.split("-")[1]))

        handle.write("@%s\n%s\n+\n%s\n" % (title, seq[:length], qual[:length]))
    handle.close()



    if(bool(RLengthTrue)):
        print("ERROR ! Dictionnary not empty.")
    else:
        print("Reads have been cut succesfully")


# ### You cannot choose the number of reads using Mitty, you can only specify the coverage of reads you want to simulate. This following script will SubSample the FastQ to have the wanted number of reads

# In[11]:


def cutting_FQ_to_wanted_reads(FileToCut,FileOut,ReadsWanted):
    
    ReadsNumber=subprocess.check_output(['wc', '-l', FileToCut])
    ReadsNumber=list(filter(None,ReadsNumber.decode("utf-8").split(" ")))
    ReadsNumber=int(ReadsNumber[0])/int(4)
    if(ReadsNumber<ReadsWanted):
        print("Error not enough reads to start with")
    else:
        NbToCut=ReadsWanted-ReadsNumber
        print(int(NbToCut))
        with open(FileToCut,'r') as file:
            lines = file.readlines()
            lines = lines[:int(NbToCut*4)]

        handle = open(FileOut, "w")
        for i in lines: 
            handle.write(i)
        handle.close()

    print("FastQ subsampled succesfully with the number of reads wanted")


# ## Simulate reads with LINEAR reference (Mitty needed)

# In[6]:


if(do_linear_simulation):
    print("Creating the reads with Mitty")
    subprocess.call("/gen_reads.sh {} {} {} {} {}".format(linear_ref_fasta,simulation_vcf_file,SampleName,coverage,Prefix), shell=True)
    print("Reads created successfully")



# ### Only keep the amount of reads that we want

# In[17]:


cutting_FQ_to_wanted_reads("{}-corrupt.fq".format(Prefix),"{}-corrupt-15Mreads.fq".format(Prefix),1500000)


# ### Now we do cut the reads into bin length size.

# In[18]:


cutting_modern_reads_to_aDNA("{}-corrupt-15Mreads.fq".format(Prefix), "{}-corrupt-15Mreads-cut.fq".format(Prefix),1500000)


# ### Transform the FastQ into Fasta to run Gargammel (Gargammel doesn't accept FastQ)

# In[19]:


subprocess.call("paste - - - - < {}-corrupt-15Mreads-cut.fq| cut -f 1,2 | sed 's/^@/>/' | tr '\t' '\n' > {}-corrupt-15Mreads-cut.fa".format(Prefix,Prefix), shell=True)
print("Transformed FastQ into Fasta successfully")


# ### Run Gargammel to add damage

# In[20]:


subprocess.call("{}/src/deamSim -mapdamage {}/examplesMapDamage/results_LaBrana/misincorporation.txt single {}-corrupt-15Mreads-cut.fa > {}-corrupt-15Mreads-cut-damaged.fa".format(PathGargammel,PathGargammel,Prefix,Prefix), shell=True)
print("Reads damaged successfully")


# ###  Add the previous baseQ (from the FastQ) to the newly created Fasta

# #### Using jupyter notebook I could run that using "!" at the beginning of the line but you probably need to run that using a bash script if you want to run it outside of jupyter notebook.

# In[15]:

get_ipython().system('cat reads.fa | sed \'s/>/@/g\' | paste - - <(seq -w 1 $(grep -c ">" reads.fa) | xargs printf \'+\\n%.s\') <(awk \'NR % 4 == 0\' reads.fq) | sed \'s/\\t/\\n/g\' > new.fq')

