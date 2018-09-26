### This program orients all sequences in a file based on the first sequence, gets rid of perfect duplicates
### and perfect partial duplicates in a fasta file and writes the remaining to a new file based on sequence length. ###

### Requres seqOrient.pl, available here http://raven.iab.alaska.edu/~ntakebay/teaching/programming/perl-scripts/seqOrient.pl
import sys
import itertools
import subprocess
import difflib
import os
#from os.path import commonprefix

def split_and_keep(s,sep):
        if not s: return [''] # consistent with string.split()
        p=chr(ord(max(s))+1)
        return s.replace(sep, p+sep).split(p)
        # http://stackoverflow.com/questions/2136556/in-python-how-do-i-split-a-string-and-keep-the-separators

class Fasta:
        def __init__(self,f):
                self.file = f
                self.filestring = f.read()
                all_lines = split_and_keep(self.filestring,">")
                all_lines = [s.split("\n",1) for s in all_lines][1:]
                self.all_lines = [[s.replace("\n","") for s in y] for y in all_lines]
                self.sorted = sorted(self.all_lines, key=lambda seq:len(seq[1]))
                self.headers = [h[0] for h in self.sorted]
                self.sequences = [s[1] for s in self.sorted]

        def isFasta(self): ## check if it's a fasta file
                if ">" == self.all_lines[0][0][0]:
                        print("Fasta file")
                else:
                        print("Not a fasta file")
                        sys.exit()

def is_dup(sequences,headers,contigID):
        uniques = []
        print(len(sequences))
        for i in range(len(sequences)-1):
                print("Sequence "+str(i)+" matches sequences "+str(i+1)+" - "+str(len(sequences)-1)+"?") 
                if sum([sequences[i] in s for s in sequences[i+1:]]) >= 1: # if sequence is in any longer ones
                        if contigID == "yes":
                                loc = [sequences[i] in s for s in sequences [i+1:]].index(True) # find first sequence it matches
                                print(headers[i].split("_")[-1]+"\t"+headers[loc].split("_")[-1])
                                print("in is_dup")
                                if headers[i].split("_")[-1] != headers[loc].split("_")[-1]: # check if contigs match
                                        print("yes, it's a duplicate, but it's on a different contig, save it")
                                        uniques.append((headers[i],sequences[i]))
                                        print(headers[i].split("_")[-1]+"\t"+headers[loc].split("_")[-1])
                                else:
                                        print("yes, it's a duplicate or partial duplicate, ignore it") # it's a duplicate, so don't save it
                                        pass
                        else:
                                print("yes, it's a duplicate or partial duplicate, ignore it") # it's a duplicate, so don't save it
                                pass
                elif sum([sequences[i] in s for s in sequences[i+1:]]) == 0: # if sequence isn't in any longer ones
                        print("no, it's unique, save it") # it's unique, so save it
                        uniques.append((headers[i],sequences[i]))
                else: # something's wrong
                        print("something's wrong with checking for duplicates at"+headers)
                        sys.exit()
 #                       print([sequences[i] in s for s in sequences[i+1:]])
        print(headers)
        print(sequences)
        print("HERE")
        if len(sequences) > 0:
                uniques.append((headers[-1],sequences[-1])) # save last sequence
        else:
                print("no sequences to append")
        return(uniques)
#       print(uniques)

def extend(uniques, contigID, k = 30): # if there's k (default = 30) characters worth of perfect overlap
        extended = []
#       k=5
        sequences = [s for h,s in uniques] # this shld be a dictionary with headers as keys
        headers = [h for h,s in uniques]
        for a,b in itertools.combinations(sequences,2):
                heada = headers[sequences.index(a)] # header a
                headb = headers[sequences.index(b)] # header b
                if len(a) > k and len(b) > k and "-" in heada and "-" in headb: # if both sequences are longer than k and neither is the reference
                        match = difflib.SequenceMatcher(None,a,b).find_longest_match(0,len(a),0,len(b))
#                       print(a)
#                       print(b)
                        print(match)
                        heada = headers[sequences.index(a)] # header a
                        headb = headers[sequences.index(b)] # header b
                        print(heada)
                        first_locations = [s.split("-") for s in heada.split("_") if "-" in s][0]
                        print(first_locations)
                        first_start = min(first_locations)
                        first_stop = max(first_locations)
                        print(headb)
                        second_locations = [s.split("-") for s in headb.split("_") if "-" in s][0]
                        print(second_locations)
                        second_start = min(second_locations)
                        second_stop = max(second_locations)
                        if second_start > first_stop and first_start < second_stop: # if they don't overlap (ie from different places on contig)
                                extended.append((headers[sequences.index(a)],a)) # save both separately
                                extended.append((headers[sequences.index(b)],b))
                                print("overlap but on different regions of contig, don't merge")
                        elif match[0] >= match[1] == 0 and match[2] > k: # if sequence a is first and overlap > k is at beginning of sequence b; or if they are perfect overlap
                                # if need contigs to match
                                # check if on same contig
                                # if on same contig, save
                                # if not on same contig, save individually
                                if contigID == "yes": # if need to be on same contig
                                        heada = headers[sequences.index(a)] # header a
                                        headb = headers[sequences.index(b)] # header b
#                                        print(heada)
#                                        print(headb)
                                        if heada.split("_")[-1] == headb.split("_")[-1]: # if on same contig
                                                merge = a[:match[0]]+b # add b to first part of sequence a
                                                extended.append((headers[sequences.index(a)]+"_"+headers[sequences.index(b)],merge)) # and save, with header of first sequence
                                                print("extended "+headers[sequences.index(a)])
                                        else:
                                                extended.append((headers[sequences.index(a)],a))
                                                extended.append((headers[sequences.index(b)],b)) 
                                                print("overlap but on different contigs, don't merge")
                                                print(heada.split("_")[-1]+"\t"+headb.split("_")[-1])
                                                pass
                                else: # if don't need to check if they're on the same contig
        #                               print(a[:match[0]])
                                        merge = a[:match[0]]+b # add b to first part of sequence a
                                        extended.append((headers[sequences.index(a)]+"_"+headers[sequences.index(b)],merge)) # and save, with header of first sequence
                                        print("extended "+headers[sequences.index(a)])
                        elif match[1] >= match[0]==0 and match[2] > k: # if sequence b is first and overlap > k is at beginning of sequence a
                                if contigID == "yes":
                                        heada = headers[sequences.index(a)] # header a
                                        headb = headers[sequences.index(b)] # header b
                                        if heada.split("_")[-1] == headb.split("_")[-1]: # if on same contig
                                                merge = b[:match[1]]+a # add a to first part of sequence b
                                                extended.append((headers[sequences.index(b)]+"_"+headers[sequences.index(a)],merge)) # and save, with header of first sequence
                                                print("extended "+headers[sequences.index(b)])
                                        else:
                                                print("overlap but on different contigs, don't merge")
                                                print(heada.split("_")[-1]+"\t"+headb.split("_")[-1])
                                                pass
                                else:
        #                               print(b[:match[0]])
                                        merge = b[:match[1]]+a # add a to first part of sequence b
                                        extended.append((headers[sequences.index(b)]+"_"+headers[sequences.index(a)],merge)) # and save, with header of first sequence
                                        print("extended "+headers[sequences.index(b)])
                        elif match[2] < k: # if no match longer than k
                                extended.append((headers[sequences.index(a)],a))
                                print("saved "+headers[sequences.index(a)])
                                extended.append((headers[sequences.index(b)],b)) # save both individually
                                print("saved "+headers[sequences.index(b)])
                        elif match[0] != 0 or match[1] != 0: # if match not at beginning of either sequence
                                extended.append((headers[sequences.index(a)],a))
                                print("saved "+headers[sequences.index(a)])
                                extended.append((headers[sequences.index(b)],b)) # save both individually
                                print("saved "+headers[sequences.index(b)])
                else: # if one or both sequences are not longer than k
                        extended.append((headers[sequences.index(a)],a))
                        print("saved "+headers[sequences.index(a)])
                        extended.append((headers[sequences.index(b)],b)) # save both individually
                        print("saved "+headers[sequences.index(b)])
#        extended_sorted = sorted(extended, key=lambda seq:seq[1])
        extended.sort(key=lambda t: len(t[1]))
        extended_sorted = extended
        print(extended_sorted)
        print([len(s[1]) for s in extended_sorted])
        return(extended_sorted)

def save_unique(uniques,var_oriented):
#       long = open("long-"+var_oriented,"a+")
#       short = open("short-"+var_oriented,"a+")
#       for i in uniques:
#               if len(i[1])<120:
#                       short.write(i[0]+"\n")
#                       short.write(i[1]+"\n")
#               elif len(i[1])>=120:
#                       long.write(i[0]+"\n")
#                       long.write(i[1]+"\n")
        filename = var_oriented.split("-")[0]+"-merged.fasta"
        save = open(filename, "a+")
        for i in uniques:
                if "reference" in i[0]:
                        pass
                else:
                        save.write(i[0]+"\n")
                        save.write(i[1]+"\n")
#                save.write(i[0]+"\n")
#                save.write(i[1]+"\n")

def main(var, contigID = "yes"):
#       var = sys.argv[1]
        var_oriented = var.split(".")[0]+"-oriented.fasta"
        if os.path.isfile(var) and os.path.getsize(var) > 0:
                print("Orienting sequences")
#                orientedfasta = subprocess.check_output(["perl", "seqOrient.pl", var], stderr=subprocess.DEVNULL) #supress output
                orientedfasta = subprocess.check_output(["perl", "seqOrient.pl", var])
                with open(var_oriented,"w") as oriented:
                        oriented.write(orientedfasta.decode("utf-8"))
                f = Fasta(open(var_oriented, "r"))
        #       f.all_lines
                f.isFasta()
        #        print(f.sorted)
        #        print(f.headers)
                print(f.sequences)
                print(var)
                uniques = is_dup(f.sequences,f.headers, contigID)
#                extended = sorted(extend(uniques, contigID), key=lambda seq:len(seq[1]))
                extended = extend(uniques, contigID)
                print(extended)
                headers = [h[0] for h in extended]
                print(len(headers))
                sequences = [s[1] for s in extended]
        #       print(len(sequences))
                new_uniques = is_dup(sequences, headers, contigID)
                save_unique(new_uniques,var_oriented)
                print("done")
        else:
                print("Can't orient sequences, file is empty or does not exist: "+var)

if __name__ == "__main__":
        main(sys.argv[1])
