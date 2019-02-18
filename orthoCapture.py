#!/bin/python3

## Test parellelization

### This program runs tblastx for a list of query sequences against a single genome,
### then BLASTs those results against the nucleotide database,
### then checks if the genes returned by the BLAST search match the identification of the query gene.
### It keeps those tblastx sequences that matched BLAST identification
### and records any genes that are missing from the initial list of query genes.
## This next part should be done by removeDups.py:
### It then checks for sequence overlap in the tblastx results that were kept and merges them. 

###Then it checks the alignment of the sequence with the query sequence
### using Water (or bl2seq?) and extends to either side if necessary. (This part has to be done by referencing [qg]_tblastx.out and blast? for that part of sequence.

# Use: python3 tblastx.py <query.csv> <genomefile.fa> <outprefix> <noskip/etc.>
### To jump to reciprocal blast: last argument should be "recblast"
### To jump to checking for match between tblastx result and BLAST identification: last argument should be "blastmatch"
### To jump to getting rid of duplicates: last argument should be "collapse"
### To jump to aligning sequences: last argument should be "align"
### To jump to extension step: last argument should be "extend"
### To run full program: last argument should be "noskip"

# The maximum reference sequence length for any given exon is 100000 nucleotides (set by pd.options.display.max_colwidth)
# Fasta headers should not have underscores, dashes

import sys
import subprocess
from subprocess import PIPE
import csv
import pandas as pd
from collections import defaultdict
import os
import itertools
import functools
import pickle
from joblib import Parallel, delayed
import multiprocessing
import removeDups

pd.options.display.max_colwidth = 100000 # set max column width so sequences aren't truncated

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return(itertools.zip_longest(*args, fillvalue=fillvalue))
        # https://docs.python.org/3.6/library/itertools.html

def csv_to_fasta(querycsv,query):
        s = csv.reader(open(querycsv, "r"))
        with open(query,"a") as o:
                for i in s:
                        if not i[1] in (None,""):
                                if i[2] in (None,""):
                                        o.write(">"+str(i[1])+"_"+str(i[3])+"\n")
                                        o.write(i[4]+"\n")
                                else:
                                        o.write(">"+str(i[1])+"_"+str(i[2])+"_"+str(i[3])+"\n")
                                        o.write(str(i[4])+"\n")

def do_add(s, x):
        l = len(s)
        s.add(x)
        return len(s) != l
        #http://stackoverflow.com/questions/27427067/python-how-to-check-if-an-item-was-added-to-a-set-without-2x-hash-lookup

def is_gene(name,genename,genevalues):
        keep = []
        count = 1
        for c in range(count):
                num_total=0
                num_matches=0
                sum_cvg=0
#                print(name+"_"+genename+'_blast.out'+" iteration "+str(c))
                for v in genevalues:
                        start = v[0]
#                        print(start)
                        stop = v[1]
#                        print(stop)
                        if os.stat(name+"/"+name+"_"+genename+'_'+start+'-'+stop+'_blast.out').st_size == 0:
                                with open(name+"/"+name+"_blast_missing.txt", "a") as m:
                                        m.write("%s\n" % genename) # everything that didn't get reciprocal blasted- everything that didn't get reciprocal blasted
                        else:
#                                print(name+'_'+genename+'_'+start+'-'+stop+'_blast.out')
                                compare = pd.read_table(name+"/"+name+'_'+genename+'_'+start+'-'+stop+'_blast.out',header=None)
                                genes = compare.iloc[:,0].str.split("_")
                                for i in range(len(genes)):
                                        gene = genes[i][0]
                                        if " "+gene+" " or " "+gene.lower()+" " or " "+gene.title()+" " or "("+gene+")" or "("+gene.lower()+")" or "("+gene.title()+")" in compare.iloc[i,2]: # if gene name is in the title of any matches
                                                sum_cvg += compare.iloc[i,1]
                                                if compare.iloc[i,1] >= 99:
                                                        num_total += 1
                                                num_matches += 1
                                        else:
                                                num_matches += 1
                                prop_total = num_total/(num_matches or 1000000) # divide by 1000000 if 0 matches
                                avg_cvg = sum_cvg/(num_matches or 1000000)
        #                        print(prop_total)
        #                        print(avg_cvg)
                                if prop_total >= .1 or avg_cvg >= 90 or num_total > 10:
                                        keep.append(name+"/"+name+"_"+genename+'_'+start+'-'+stop+'_blast.out')
                                        # keep longest and best matches
        with open(name+"/"+name+"_tblastx_kept.out", "a") as k:
                for item in keep:
                        k.write("%s\n" % item)

def is_missing(name,querycsv):
        allgenes = pd.read_csv(querycsv,dtype=object)
#        querynames = [el[1] for el in allgenes.iloc[:,0].str.split(" ")] # in query
        allgenenames = allgenes.iloc[:,1].astype(str).str.cat(allgenes.iloc[:,2].astype(str).str.cat(allgenes.iloc[:,3].astype(str), sep="_"), sep="_")
#        print(allgenenames.head())
        querynames = set(allgenenames)
        missing = []
        with open(name+"/"+name+"_tblastx_kept.out","r") as l:
                keptlines = l.read().splitlines()
        for n in querynames:
                if any(n not in k for k in keptlines):
                        missing.append(n)
#        print(missing)
        with open(name+"/"+name+"_tblastx_missing.out", "a") as m:
                for item in missing:
                        m.write("%s\n" % item)

def combine_kept(name,genename,genevalues): # makes fasta file for each gene that has a sequence in this reference genome
#        strip = "_"
#        geneonly = strip.join(genename.split(strip)[:-1])
        with open(name+"/"+name+"_tblastx_kept.out", "r") as k:
                filenames = k.read().splitlines()
        with open(name+"/"+name+"_"+genename+"_match_all_tblastx.fa", "a") as outfile:
                print("Creating "+name+"/"+name+"_"+genename+"_match_all_tblastx.fa")
                for i in genevalues:
                        genefiles = [s.rsplit("_",1)[0]+"_tblastx.fa" for s in filenames if genename in s]
                        for fname in genefiles:
#                                print(fname)
#                                print(name)
                                with open(fname) as infile:
                                        for line in infile:
                                                outfile.write(line)

def add_ref(name,querycsv):
#        print(name)
        allquery = pd.read_csv(querycsv,dtype=object)
#        print(allquery.head())
        with open(name+"/"+name+"_tblastx_kept.out","r") as f:
                rows = f.read().splitlines()
        already = []
        for row in rows:
                if len(row.split("_")[:-2]) == 3:
                        genename, ref = row.split("_")[1:-2]
                        exon = "nan"
#                        print(exon)
                elif len(row.split("_")[:-2]) == 4:
                        genename, exon, ref = row.split("_")[1:-2]
#                        print(exon+" exon for add_ref")
                else:
                        print("There's an issue in add_ref.")
                        print(row.split("_")[:-2])
                        sys.exit()
                if name and genename and ref:
#                        print(name)
                        if exon == "nan":
                                filename = name+"/"+name+"_"+genename+"_"+ref+"_match_all_tblastx.fa"
                        else:
                                filename = name+"/"+name+"_"+genename+"_"+exon+"_"+ref+"_match_all_tblastx.fa"
                        if filename in already:
                                pass
#                                print("Already put a reference sequence in "+filename)
                        elif filename not in already:
                                try:
                                        seq = allquery[(allquery["Gene Name Shorthand"] == genename) & (allquery["Exon Number"].apply(str) == str(exon)) & (allquery["Design Species"].str.contains(ref))]["Design Sequence"].to_string().split()[1] # convert pandas entry to string, get rid of row identifier
#                                        print(seq)
#                                        print(filename)
                                        with open(filename, "r") as infile:
                                                alltext = infile.read()
                                        with open(filename, "w") as outfile:
                                                outfile.write(">"+ref+" reference\n")
                                                outfile.write(seq+"\n"+alltext)
                                        already.append(filename)
                                        with open(name+"/"+name+"_amended.out", "a") as am:
                                                am.write(filename+"\n")
#                                        print(filename)
#                                        print("File amended")
                                except IOError:
                                        print("Uh oh, can't open file "+filename)
                        else:
                                print("Something's wrong at "+filename)
                else:
                        print("Something's missing: "+name+" "+genename+" "+ref)
        return(already)

def align(name,filename,querycsv):
        genename = "_".join(filename.split("_")[1:-3])
        geneshort = genename.split("_")[0]
        reference = genename.split("_")[-1]
        if genename.split("_")[-1] == genename.split("_")[1]:
                exonnumber = "nan"
        else:
                exonnumber = genename.split("_")[1]
        allquery = pd.read_csv(querycsv,dtype=object)
        subject = allquery[(allquery["Gene Name Shorthand"] == geneshort) & (allquery["Exon Number"].apply(str) == str(exonnumber)) & (allquery["Design Species"].str.contains(reference))]["Design Sequence"].to_string().split()[1]
        with open(genename+"_reference.fasta","w") as ref:
                ref.write(">"+reference+"\n")
                ref.write(subject+"\n")
        if os.path.exists(name+"/"+name+"_"+genename+"_match_all_tblastx-merged.fasta"):
                with open(name+"/"+name+"_"+genename+"_match_all_tblastx-merged.fasta", "r") as allfile:
                        print("aligning "+name+"/"+name+"_"+genename+"_match_all_tblastx-merged.fasta")
                        lines = allfile.readlines()
                        pairs = grouper(lines,2)
                        start=[]
                        stop=[]
                        length=[]
                        header=[]
                        sequence=[]
                        for i in pairs:
                                a,b = i
                                header.append(a)
                                sequence.append(b)
                                if a[0] != ">":
                                        print("Not a proper fasta file at "+a)
                                        sys.exit()
                                else:
                                        locations = [s.split(" ")[0] for s in a.split("_") if "-" in s] # look for list elements with a dash since that's the start-stop separator, split location from "Reverse:" if present
                                                # this would be a problem if the sequences were in two different orientations, but if they're on the same contig they shouldn't be
                                        if len(locations) == 2: # if it was merged st header has two sequence names
#                                                print(str(len(locations))+" locations")
                                                first_start = min(min([locations[0].split("-")]))
                                                second_start = min(min([locations[1].split("-")]))
                                                first_stop = max(max([locations[0].split("-")]))
                                                second_stop = max(max([locations[1].split("-")]))
                                                if second_start > first_stop or first_start > second_stop: # if they don't overlap (ie from different places on contig)
#                                                if second_start > first_stop and first_start < second_stop: # if they don't overlap (ie from different places on contig)
                                                        start.append(first_start) # add them separately
                                                        stop.append(first_stop)
                                                        start.append(second_start)
                                                        stop.append(second_stop)
                                                else:
#                                                        #str(len(locations))+" locations")
                                                        start.append(min(min([pos.split("-") for pos in locations]))) # smallest location value
                                                        stop.append(max(max([pos.split("-") for pos in locations]))) # largest location value                                       
                                        elif len(locations) == 1:
#                                                print(str(len(locations))+" locations")
                                                start.append(min(min([pos.split("-") for pos in locations]))) # smallest location value
                                                stop.append(max(max([pos.split("-") for pos in locations]))) # largest location value                                                        stop.append(max(max([pos.split("-") for pos in locations]))) # largest location value
                        for j in range(len(start)): # this aligns only the longer section
                                length.append(int(stop[j])-int(start[j]))
#                                print(stop[j]+" - "+start[j])
#                        print(length)
                        best = length.index(max(length))
                        if max(length) < 120:
                                with open(name+"/"+name+"_tblastx.short","a") as log:
                                        log.write("Careful, the longest alignment for "+genename+" is only "+str(max(length))+" nucleotides.\n")
                        bestheader = header[best]
                        bestsequence = sequence[best]
                        headersplit = bestheader.split(">")[1:] # split on header carat (and get rid of resultant first empty entry), will have one entries if one carat, two or more if multiple carats (which indicates multiple sequences merged)
                        separate = [s.split("_") for s in headersplit] # separate headers
                        contigs = [[p[p.index(a)+1].strip() for a in p if "-" in a] for p in separate] # get element after start-stop (since that's the uniquely identifiable one)
#                        print(contigs[0][0].split()[0]) # have to split contigs[0][0] to get rid of " Reversed:" if exists
                        samecontigs = all(x == contigs[0] for x in contigs)
                        if samecontigs: # if all contigs the same
                                with open(name+"/"+name+"_"+genename+"_temp.fa", "w") as temp:
                                        temp.write(bestheader)
                                        temp.write(bestsequence+"\n")
                                tblastx_command_line = 'ncbi-blast-2.2.29+/bin/tblastx -query '+genename+"_reference.fasta"+' -subject '+name+"/"+name+"_"+genename+"_temp.fa"+' -out '+name+"/"+name+'_'+genename+"_"+start[best]+"-"+stop[best]+"_"+contigs[0][0].split()[0]+'_align.out -outfmt "6 sseqid qstart qend sstart send qcovhsp slen qlen bitscore evalue" -max_target_seqs 1'
                                        # using ncbi blast version 2.2.29 so qcovhsp outputs percent nucleotide coverage
                                        # coverage is amount of reference sequence aligned to pulled sequence
                                subprocess.Popen(tblastx_command_line, shell=True).wait() # shell because of argument in quotes, wait so script doesn't advance until done
                                        # will output file [name]_align.out
                                with open(name+"/"+name+"_aligned.out", "a") as al:
                                        al.write(name+"/"+name+'_'+genename+"_"+start[best]+"-"+stop[best]+"_"+contigs[0][0].split()[0]+'_align.out\n')
                        else:
                                print("merged sequences are from different contigs, need to extend separately on each side")
        else:
                print("Merged file does not exist, cannot align: "+filename)

def extend(name, filename, querycsv, genome):
        locations = [s for s in filename.split("_") if "-" in s][0].split("-")
        start = int(locations[0])
        stop = int(locations[1])
        if os.path.exists(filename) and os.path.getsize(filename) > 0:
                alignments = pd.read_table(filename,header=None)
                qlen = alignments[alignments[9] == alignments[9].min()].iloc[0,7]
                slen = alignments[alignments[9] == alignments[9].min()].iloc[0,6]
                sstart = min([alignments[alignments[9] == alignments[9].min()].iloc[0,3],alignments[alignments[9] == alignments[9].min()].iloc[0,4]])
                send = max([alignments[alignments[9] == alignments[9].min()].iloc[0,3],alignments[alignments[9] == alignments[9].min()].iloc[0,4]])
                qstart = min([alignments[alignments[9] == alignments[9].min()].iloc[0,1],alignments[alignments[9] == alignments[9].min()].iloc[0,2]])
                qend = max([alignments[alignments[9] == alignments[9].min()].iloc[0,1],alignments[alignments[9] == alignments[9].min()].iloc[0,2]])
                cvg = slen/qlen * 100
        #        first=alignments[alignments[9] == alignments[9].min()].iloc[0,3] # best evalue is best alignment
        #        last=alignments[alignments[9] == alignments[9].min()].iloc[0,4]
                if cvg < 100:
                        sstartpos = [alignments[alignments[9] == alignments[9].min()].iloc[0,3],alignments[alignments[9] == alignments[9].min()].iloc[0,4]].index(sstart)
                        qstartpos = [alignments[alignments[9] == alignments[9].min()].iloc[0,1],alignments[alignments[9] == alignments[9].min()].iloc[0,2]].index(qstart)
                        if sstartpos == 0 and qstartpos == 0: # both running forward
                                sfrontdiff = sstart - 1 # difference between subject s$
                                        # if > 0, need to subtract from difference bet$
                                senddiff = slen - send # difference between subject le$
                                        # if > 0, need to subtract from difference bet$
                                qfrontdiff = qstart - 1 # difference between query sta$
                                qenddiff = qlen - qend # difference between query leng$
                                frontmissing = qfrontdiff - sfrontdiff # since running$
                                endmissing = qenddiff - senddiff # since running forwa$
                                newstart = max((start - frontmissing),1)
                                newstop = stop + endmissing
                        elif sstartpos == 0 and qstartpos == 1: # subject running forward, query (reference) running backwa$
                                sfrontdiff = sstart - 1
                                senddiff = slen - send
                                qfrontdiff = qlen - qend
                                qenddiff = qstart - 1
                                frontmissing = qfrontdiff - sfrontdiff
                                endmissing = qenddiff - senddiff
                                newstart = max((start - frontmissing),0)
                                newstop = stop + endmissing
                        elif sstartpos ==1 and qstartpos == 0: # subject running backward, query (reference) running forward
                                sfrontdiff = slen - send
                                senddiff = sstart - 1
                                qfrontdiff = qstart - 1
                                qenddiff = qlen - qend
                                frontmissing = qfrontdiff - sfrontdiff
                                endmissing = qenddiff - senddiff
                                newstart = max((start - endmissing),0)
                                newstop = stop + frontmissing
                        elif sstartpos == 1 and qstartpos == 1: # both running backward
                                sfrontdiff = sstart - 1 # difference between subject s$
                                        # if > 0, need to subtract from difference bet$
                                senddiff = slen - send # difference between subject le$
                                        # if > 0, need to subtract from difference bet$
                                qfrontdiff = qstart - 1 # difference between query sta$
                                qenddiff = qlen - qend # difference between query leng$
                                frontmissing = qfrontdiff - sfrontdiff # since running$
                                endmissing = qenddiff - senddiff # since running backw$
                                newstart = max((start - frontmissing),1)
                                newstop = stop + endmissing
                        header = alignments[alignments[9] == alignments[9].min()].iloc[:,0].to_string().split()[1]
        #                print(header)
                        headersplit = header.split(">")[0:] # split on header carat (and get rid of resultant first empty entry), will have one entries if one carat, two or more if multiple carats (which indicates multiple sequences merged)
        #                print(header.split(">"))
                        separate = headersplit[0].split("_")
        #                print(headersplit)
        #                print(separate)
        #[s.split("_") for s in headersplit] # separate headers
                        contig = [separate[separate.index(a)+1] for a in separate if "-" in a][0]
        #                print(contig)
        #                print(type(newstart))
        #                print(type(newstop))
                        blastdb_command_line = 'ncbi-blast-2.7.1+/bin/blastdbcmd -db '+genome+' -entry '+contig+' -range '+str(newstart)+'-'+str(newstop)
        #                print(blastdb_command_line)
                        out = subprocess.Popen(blastdb_command_line, shell=True, stdout=PIPE) # to output so can modify header
                        output = out.stdout.read()
                        sequence = "".join([s.decode("utf-8") for s in output.splitlines()[1:]])
#                        print(filename.rsplit("_",3)[:-3][0])
                        headerextend = ">"+filename.rsplit("_",3)[:-3][0]
        #                print(sequence)
#                        print(headerextend)
                        with open(name+"/"+name+"_extended.fa","a") as f:
                                f.write(headerextend+"_"+str(newstart)+"-"+str(newstop)+"_"+contig+"\n")
                                f.write(sequence+"\n")
                elif cvg == 100:
                        fullstart = min(start,stop)
                        fullstop = max(start,stop)
                        header = alignments[alignments[9] == alignments[9].min()].iloc[:,0].to_string().split()[1]
#                        print(header)
                        headersplit = header.split(">")[0:] # split on header carat (and get rid of resultant first empty entry), will have one entries if one carat, two or more if multiple carats (which indicates multiple sequences merged)
                        separate = headersplit[0].split("_")
                        contig = [separate[separate.index(a)+1] for a in separate if "-" in a][0]
                        blastdb_command_line = 'ncbi-blast-2.7.1+/bin/blastdbcmd -db '+genome+' -entry '+contig+' -range '+str(fullstart)+'-'+str(fullstop)
                        out = subprocess.Popen(blastdb_command_line, shell=True, stdout=PIPE) # to output so can modify header
                        output = out.stdout.read()
                        sequence = "".join([s.decode("utf-8") for s in output.splitlines()[1:]])
                        with open(name+"/"+name+"_extended.fa", "a") as f:
                                f.write(">"+name+"_"+header+"\n")
                                f.write(sequence+"\n")
#                        print("full coverage without extending")
        else:
                print("Alignment does not exist, cannot extend: "+filename)

def tblastx(name,querycsv,query,genome):
        if not os.path.isfile(query):
               csv_to_fasta(querycsv,query)
        tblastx_command_line = 'ncbi-blast-2.2.29+/bin/tblastx -query '+query+' -db '+genome+' -culling_limit 1 -out '+name+'/'+name+'_tblastx.out -outfmt "6 qseqid qstart qend sstart send qcovhsp bitscore sacc" -max_target_seqs 1'
                # using ncbi blast version 2.2.29 so sseq outputs nucleotide sequence (2.6.* outputs protein sequence)
        print("tblastx command line")
        print(tblastx_command_line)
        subprocess.Popen(tblastx_command_line, shell=True).wait() # shell because of argument in quotes, wait so script doesn't advance until done
                # will output file [name]_tblastx.out
        f = open(name+'/'+name+"_tblastx.out","r")
        reader = csv.reader(f, delimiter = "\t")
        fasta = defaultdict(list)
        for i in reader:
                fasta[i[0]].append([i[3],i[4],i[7]]) # gene into dictionary with start, stop, sequence, and contig as values
        f.close()
        with open(name+'/'+name+'_fasta_tblastx.pickle', 'wb') as handle:
                pickle.dump(fasta, handle, protocol=pickle.HIGHEST_PROTOCOL)

def recBlast(fasta,gn,name,genome):
        if ntloc="NONE":
                sys.exit("You didn't specify the location of your nucleotide database. Run prep_directory.sh <ntdb location> and try again from the recblast step.")
        elif not os.path.isdir(ntloc):
                sys.exit("The file you specified for the location of your nucleotide database doesn't exist. Check the path and run ntdb_location.sh <ntdb location> and try again from the recblast step.")
        else:
                for x in range(len(fasta[gn])): # for each instance of that gene
                        filename = name+"/"+name+"_"+gn+"_"+fasta[gn][x][0]+'-'+fasta[gn][x][1]+"_tblastx.fa"
                        with open(filename,"w") as n: # create a new file with the name and sequence
                                start = min(fasta[gn][x][0], fasta[gn][x][1])
                                stop = max(fasta[gn][x][0], fasta[gn][x][1]) # because blastdbcmd can't handle start greater than stop
                                n.write(">"+gn+"_"+start+"-"+stop+"_"+fasta[gn][x][2]+"\n")
                                blastdb_command_line = 'ncbi-blast-2.7.1+/bin/blastdbcmd -db '+genome+' -entry '+fasta[gn][x][2]+' -range '+str(start)+'-'+str(stop)
                                print("blastdbcmd command line:")
                                print(blastdb_command_line)
                                out = subprocess.Popen(blastdb_command_line, shell=True, stdout=PIPE)
                                output = out.stdout.read()
                                sequence = "".join([s.decode("utf-8") for s in output.splitlines()[1:]]) # decode wrapped sequence as list elements from bytes to string, then join list elements
#                                print(sequence)
                                fasta[gn][x].append(sequence) # add sequence
                                n.write(sequence+"\n")
#                        print(filename)
                        blast_command_line = 'ncbi-blast-2.7.1+/bin/blastn -task dc-megablast -db '+ntloc+' -query '+filename+' -max_hsps 50 -out '+name+"/"+name+"_"+gn+"_"+fasta[gn][x][0]+'-'+fasta[gn][x][1]+'_blast.out -outfmt "6 qseqid qcovhsp stitle"'
                        subprocess.Popen(blast_command_line, shell=True).wait() # blast that sequence
                        print("blastn command line:")
                        print(blast_command_line)

def par_recBlast(name,genome):
        num_cores = 20
#        num_cores = multiprocessing.cpu_count()
        with open(name+'/'+name+'_fasta_tblastx.pickle', 'rb') as handle:
                fasta = pickle.load(handle)
        Parallel(n_jobs=num_cores)(delayed(recBlast)(fasta,gn,name,genome) for gn in fasta.keys())
        with open(name+'/'+name+'_fasta_recblast.pickle', 'wb') as handle:
                pickle.dump(fasta, handle, protocol=pickle.HIGHEST_PROTOCOL)

def par_blastMatch(name,querycsv):
#        num_cores = multiprocessing.cpu_count()
        num_cores = 20
        with open(name+'/'+name+'_fasta_recblast.pickle', 'rb') as handle:
                fasta = pickle.load(handle)
        for i,j in fasta.items():
#                newi="_".join(i.split("_")[:-1])
                is_gene(name,i,j)
        is_missing(name,querycsv)
        for i,j in fasta.items():
#                newi="_".join(i.split("_")[:-1])
                combine_kept(name,i,j)
        add_ref(name,querycsv)
#        for i in fasta.keys():
#                combine_kept(name,i,len(fasta[i]))
#                removeDups.main(name+"/"+name+"_"+i+"_match_all_tblastx.fa","yes")
        Parallel(n_jobs=num_cores)(delayed(removeDups.main)(name+"/"+name+"_"+i+"_match_all_tblastx.fa","yes") for i in fasta.keys())

def collapse(name,querycsv):
        with open(name+'/'+name+'_fasta_recblast.pickle', 'rb') as handle:
                fasta = pickle.load(handle)
        add_ref(name,querycsv)
        for i,j in fasta.items():
#                print("collapse "+i)
                removeDups.main(name+"/"+name+"_"+i+"_match_all_tblastx.fa","yes")

def main():
#        print(sys.version)
        querycsv = sys.argv[1]
        genome = sys.argv[2]
        name = sys.argv[3]
        query = querycsv.split(".")[0]+".fa"
        skip = sys.argv[4]
        # Define ntloc
        ntloc = '/test/test/path'
        if skip == "noskip":
                if not os.path.isdir(name):
                        os.makedirs(name) # make a directory for all of the files
                else:
                        print(name+" directory already exists")
                tblastx(name,querycsv,query,genome)
                par_recBlast(name,genome)
                par_blastMatch(name,querycsv)
                collapse(name,querycsv)
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                for i in already:
                        if os.path.exists(i):
                                align(name, i, querycsv)
                with open(name+"/"+name+"_aligned.out","r") as al:
                        aligned=al.read().splitlines()
                for i in aligned:
                        extend(name, i, querycsv, genome)
        if skip == "tblastx":
                if not os.path.isdir(name):
                        os.makedirs(name) # make a directory for all of the files
                else:
                        print(name+" directory already exists")
                tblastx(name,querycsv,query,genome)
                print("tblastx step complete. go to recblast.")
                sys.exit("tblastx step complete. go to recblast.")
        if skip == "recblast":
                par_recBlast(name,genome)
                print("recblast step complete. go to blastmatch.")
                sys.exit("recblast step complete. go to blastmatch.")
                par_blastMatch(name,querycsv)
                collapse(name,querycsv)
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                for i in already:
                        if os.path.exists(i):
                                align(name, i, querycsv)
                with open(name+"/"+name+"_aligned.out","r") as al:
                        aligned=al.read().splitlines()
                for i in aligned:
                        extend(name, i, querycsv, genome)
        if skip == "blastmatch":
                par_blastMatch(name,querycsv)
                print("blastmatch step complete. go to collapse.")
                sys.exit("blastmatch step complete. go to collapse.")
                collapse(name,querycsv)
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                for i in already:
                        if os.path.exists(i):
                                align(name, i, querycsv)
                with open(name+"/"+name+"_aligned.out","r") as al:
                        aligned=al.read().splitlines()
                for i in aligned:
                        extend(name, i, querycsv, genome)
                sys.exit()
        if skip == "collapse":
                collapse(name,querycsv)
                print("collapse step complete. go to align.")
                sys.exit("collapse step complete. go to align.")
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                for i in already:
                        if os.path.exists(i):
                                align(name, i, querycsv)
                with open(name+"/"+name+"_aligned.out","r") as al:
                        aligned=al.read().splitlines()
                for i in aligned:
                        extend(name, i, querycsv, genome)
        if skip == "align":
                with open(name+"/"+name+"_amended.out","r") as am:
                        already=am.read().splitlines()
                for i in already:
                        if os.path.exists(i):
                                align(name, i, querycsv)
                print("align step complete. go to extend.")
                sys.exit("align step complete. go to extend.")
                with open(name+"/"+name+"_aligned.out","r") as al:
                        aligned=al.read().splitlines()
                for i in aligned:
                        extend(name, i, querycsv, genome)
        if skip == "extend":
                with open(name+"/"+name+"_aligned.out","r") as al:
                        aligned=al.read().splitlines()
                for i in aligned:
                        extend(name, i, querycsv, genome)
                print("extend step complete. done.")

if __name__ == "__main__":
        main()
