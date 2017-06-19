#!/usr/bin/python

'''
Go through Illlumina MEGA manifest file and transform mixed classes of identifiers to unified Chr:Pos-Ref-Alt format
(simplifies merging with ExAC, ClinVar, 1000 Genomes etc) 

 *  replacing monomorphic SNP IDs that have a "0" allele in the MEGA IDs with the correct allele
 *  replacing indels that have a D or I allele in the MEGA IDs with the correct base(s)
 *  avoid PLUS and MINUS confusion for indels by using all plus-oriented (to match ClinVar)
 *  avoid TOP and BOT confusion by making all top-oriented (to match ClinVar) --> May need to reexamine these later.

 *  if a ref-alt switch of an INDEL or SNP matches to Clinvar: Switch them, but flag it
 *  if a reverse complement of a SNP matches to Clinvar: take the rev comp, but flag it

Takes in the Illumina MEGA manifest without headers, the Clinvar file, and a list of unique clinvar IDs in format chr: pos- ref -alt. 

Returns a CSV file with the old MEGA ID and new (Clinvar-matching) MEGA ID, and other columns as needed. 


'''

from Bio import SeqIO

fastafile = '/path/to/MEGA_Consortium_15063755_B2.fa' #Made using bedtools getfasta on 4-25-16

fasta_sequences = SeqIO.parse(open(fastafile), 'fasta')

refids ={  ":".join( [fasta.id.split(":")[0], fasta.id.split("-")[1] ]): str(fasta.seq) for fasta in fasta_sequences } 

#csvfile= 'temp.csv'
csvfile = '/path/to/MEGA_Consortium_15063755_B2.noheader.csv'

outfilename = 'result.txt'

outfile = file(outfilename, 'w')

errfilename= 'err.txt'
errfile = file(errfilename,'w')

clinvar = '/path/to/clinvar.tsv'

clinvar_id_list = '/path/to/clinvar.list'

clinvar_positions = set( [':'.join([line.split()[0], line.split()[1]])  for line in file(clinvar) ])

rc = {"A":"T", "G":"C", "T":"A","C":"G"}

clinvar_ids = set([line.strip() for line in file(clinvar_id_list)])

counter = 0


for line in file(csvfile):
    name,chrom,pos,allele,strand,ref,alt,cat,refid, real_ref,flag='','','','','','','','','','','flag0'
    line=line.strip().split(',')
    name = line[1]
    chrom,pos = line[9:11]
    allele = line[3]
    strand = line[2] # For indels, use the Illumina strand 
    flag = 'flag0' 
    if allele.split("[")[1].startswith("D"): # If it's an indel 
        cat= "indel"
        # If there are bases in the Name (Col2):  use them -- the D-I order is often wrong 
        if name.endswith( ('A', 'C', 'T','G')  ):
            ref, alt = name.split("-")[1], name.split("-")[2] 
  
        else: # Concat the last letter before bracket and all the stuff after - to make ref allele
            ref = line[17].split('[')[0][-1] + line[17].split('/')[1].split("]")[0] 
            alt = line[17].split('[')[0][-1]
            if strand=="PLUS":
                pass
                
            elif strand == "MINUS":
                # Switch ref and alt
                ref, alt = alt, ref 


        if '-'.join([':'.join([chrom,pos]), ref, alt]) not in clinvar_ids:
            
            # If switching ref and alt alleles finds a match 
            if '-'.join([':'.join([chrom,pos]), alt, ref]) in clinvar_ids: 
                #print 'Swithcng these alleles:\n%s' % ('-'.join([':'.join([chrom,pos]), ref, alt]))
                alt, ref = ref, alt
                flag ='flag1'

    elif allele.split("[")[1].startswith("I"): #If labelled as insertion
        cat= "indel"
           # If there are bases in the Name (Col2) use them because the D-I order is often wrong
        if name.endswith( ('A', 'C', 'T','G')  ):
            ref, alt = name.split("-")[1], name.split("-")[2] 


        else: # Concat the last letter before bracket and all the stuff after - to make ref allele
            ref = line[17].split('[')[0][-1]
            alt = line[17].split('[')[0][-1] + line[17].split('/')[1].split("]")[0]
                
            if strand == "MINUS":
                # Switch ref and alt
                ref, alt = alt, ref 

        if '-'.join([':'.join([chrom,pos]), ref, alt]) not in clinvar_ids:
            
            # See if the alleles are switched 
            if '-'.join([':'.join([chrom,pos]), alt, ref]) in clinvar_ids: 
                #print 'Switching these alleles:\n%s' % ('-'.join([':'.join([chrom,pos]), ref, alt]))
                alt, ref = ref, alt
                flag ='flag1'

    # If a SNP
    elif allele.split("[")[1].startswith("A") or allele.split("[")[1].startswith("C") or allele.split("[")[1].startswith("T") or allele.split("[")[1].startswith("G"):
    #elif ':'.join ([chrom, pos]) == "17:62018773": #TODO: Case where RSIDs not correct in MEGA but can be used to get correct allele
    #elif ":".join([chrom,pos]) == "X:79947587":
        #break
    #elif ":".join ([chrom, pos])  == "12:13906449": #Example of case that should not overlap with MEGA bc its a diff ClinVar allele
        cat = "snp"
        strand = line[15] #Source strand info more useful in this case 
        ref = line[17].split("]")[0][-1]
        alt = line[17].split("[")[1][0]
        
        #### Get the ref allele from refids
        try:
            real_ref = refids[':'.join([chrom,pos] )]
            if ref == real_ref: # If hg19 base is 1st base in brackets, ref: do nothing
                pass
            elif alt == real_ref: 
                flag='flag1'
                ref, alt = alt, ref  # If hg19 is 2nd base in brackets, alt: need to switch ref and alt
            elif rc[ref] == real_ref: 
                flag='flag2'
                ref, alt = rc[ref], rc[alt]
            elif rc[alt] == real_ref:
                flag='flag3'
                ref, alt = rc[alt], rc[ref]

        except KeyError: # Chrom 0 , pos 0 [can be rescued for ClinGen variants--  see below, flag 4]
            cat = "no-pos"
            flag = 'flag-1'
 
        # if strand == "BOT":
        #     # Taka the reverse complement 
        #     ref, alt = rc[ref], rc[alt]
        
        else:
            pass

        if '-'.join([':'.join([chrom,pos]), ref, alt]) in clinvar_ids:
            flag = 'flag0'
        # If this allele is not in ClinVar 
        else:
            
            # See if the rev comp is in clinvar
            if '-'.join([':'.join([chrom,pos]), rc[ref], rc[alt]]) in clinvar_ids:
                #print 'Taking rev comp for %s' % ('-'.join([':'.join([chrom,pos]), rc[ref], rc[alt]]))
                ref = rc[ref]
                alt = rc[alt]
                flag ='flag4'   
            
            #elif '-'.join([':'.join([chrom,pos]), alt, ref]) in clinvar_ids: 
                #print 'Switching these alleles:\n%s' % ('-'.join([':'.join([chrom,pos]), ref, alt]))
                #alt, ref = ref, alt
                #flag = 3 #FIXME: I think this one is dangerous. It may result in mistakes! 
    else:
        #print 'error: %s is neither indel nor SNP   ' % (line)
        cat = "unk"

    # writing out

    # If it's in Clinvar positions, write it, else skip for now. 

    if ':'.join( [chrom, pos]) in clinvar_positions:
        
        # If the positions match but this entry is not in the Clinvar IDs

        if '-'.join([':'.join([chrom,pos]), ref, alt]) not in clinvar_ids:
            # If the original MEGA name is in ClinVar IDs (a small number of cases)
            
            if name in clinvar_ids:

                chrom,pos,ref, alt = name.split(":")[0], name.split(":")[1].split("-")[0], name.split(":")[1].split("-")[1], name.split(":")[1].split("-")[2]
                flag = 'flag5'
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), flag,cat)) 

            else:
                flag= 'flag9'
                outfile.write ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), flag,cat))
            # Look at whether the original ID is actually in the 
        elif '-'.join([':'.join([chrom,pos]), ref, alt]) in clinvar_ids:
            # Check that it's the same ref in refseq
            if ref == refids[':'.join([chrom,pos] )]:

                #outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), ','.join([line[2], line[3], line[13], line[15], line[16], line[17] ]) ))
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), flag, cat)) 
            else: 
                pass
        else:
            pass
    else: # If it's not in Clinvar IDs, let's fix it anyway. 
        if cat == "snp":

            if ref == refids[':'.join([chrom,pos] )]: #If it's in hg19, and same strand, ref allele is already correct
                flag='flag0'
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), flag,cat)) 
            elif rc[ref]== refids[':'.join([chrom,pos]) ]: # If it's in hg19, and reverse strand, take the rev comp of this position
                flag='flag6'
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, rc[ref], rc[alt], '-'.join([':'.join([chrom,pos]), rc[ref], rc[alt]]), flag,cat)) 
            else: 
                flag = 'flag7'
                outfile.write ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]),  flag, cat ))

        if cat == "no-pos": # Try to use the name to help out 
            try:
                chrom,pos,ref, alt = name.split(":")[0], name.split(":")[1].split("-")[0], name.split(":")[1].split("-")[1], name.split(":")[1].split("-")[2]
                flag = 'flag-1'
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), flag,cat)) 
            except IndexError:
                errfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), ','.join([line[2], line[3], line[13], line[15], line[16], line[17] ]) ))
        elif cat == "indel": #if it's an indel that wasnt in ClinVar
            flag = 'flag8'
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, chrom, pos, ref, alt, '-'.join([':'.join([chrom,pos]), ref, alt]), flag,cat)) 
        else: 
            pass
outfile.close()

