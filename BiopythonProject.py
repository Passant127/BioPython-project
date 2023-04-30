import os
import sys
import getopt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# --------------------------------------------------------------------

def gc(dna):
    nbases = dna.count('N')+dna.count('n')
    gcpercent = float(dna.count('G')+dna.count('g') +
                      dna.count('C')+dna.count('c'))/(len(dna)-nbases)
    return gcpercent
# --------------------------------------------------------------------


def is_valid(seq, typee):
    flag = 1
    typee = typee.lower()
    seqUpper = seq.upper()
    seqUpperList = list(seqUpper)
    seqLen = len(seqUpperList)
    valid_DNA = ['A', 'G', 'C', 'T']
    valid_RNA = ['A', 'G', 'C', 'U']
    valid_Protein = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                     'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    if typee == 'protein':
        for i in range(seqLen):
            if seqUpperList[i] in valid_Protein:
                pass
            else:
                flag = 0

    elif typee == 'dna':
        for i in range(seqLen):
            if seqUpperList[i] in valid_DNA:
                pass
            else:
                flag = 0

    elif typee == 'rna':
        for i in range(seqLen):
            if seqUpperList[i] in valid_RNA:
                pass
            else:
                flag = 0

    if flag == 0:
        return False
    else:
        return True
# --------------------------------------------------------------------


def transcribe(dna):
    if (is_valid(dna, 'dna') == True):
        dna = dna.upper()
        for i in range(len(dna)):
            if dna[i] == 'T':
                x = dna.replace(dna[i], 'U')
        print(x)
    else:
        print("NOT VALID SEQUENCE")
# --------------------------------------------------------------------


def seq_alignment(seq1, seq2, *file):
    try:
        if (is_valid(seq1, 'dna') == True and is_valid(seq2, 'dna') == True):
            f = open(file[0], 'a')
            alignments = pairwise2.align.globalxx(seq1, seq2)
            for alignment in alignments:
                f.write(format_alignment(*alignment))
        else:
            print("NOT VALID SEQUENCES")
    except:
        if (is_valid(seq1, 'dna') == True and is_valid(seq2, 'dna') == True):
            alignments = pairwise2.align.globalxx(seq1, seq2)
            for alignment in alignments:
                print(format_alignment(*alignment))
        else:
            print("NOT VALID SEQUENCES")
# --------------------------------------------------------------------


def seq_alignment_files(file1, file2, *file3):
    try:
        seqoffile1 = SeqIO.read(file1, "fasta")
        seqoffile2 = SeqIO.read(file2, "fasta")
    except:
        print("An error occurred please check that you entered to files of type fasta file with a valid names and contain one sequence")
    try:
        align = pairwise2.align.globalxx(
            seqoffile1, seqoffile2)  # to align to sequences
        f = open(file3[0], 'w')
        for alignment in align:
            # changed to string to be write on file
            alignmentt = str(alignment)
            f.write(alignmentt)  # not in a format
            formatalign = str(format_alignment(*alignment))  # in a format
            f.write(formatalign)
            f.close()
    except:
        align = pairwise2.align.globalxx(seqoffile1, seqoffile2)
        for alignment in align:
            print(alignment)
            print(format_alignment(*alignment))

# --------------------------------------------------------------------


def merge_fasta(file1, file2, output="", *files):
    try:
        seqoffile1 = SeqIO.read(file1, "fasta")
        seqoffile2 = SeqIO.read(file2, "fasta")
    except:
        print("An error occurred please check that you entered to files of type fasta file with a valid names")
    Seqs = []
    Seqs.append(seqoffile1)
    Seqs.append(seqoffile2)
    if len(files) != 0:
        for file in files:
            try:
                z = SeqIO.read(file, "fasta")
                Seqs.append(z)
            except:
                print(
                    "An error occurred please check that you entered to files of type fasta file with a valid names")
    if output != "":
        with open(output, "w") as fw:
            for i in Seqs:
                SeqIO.write(i, fw, 'fasta')
            fw.close()
    if output == "":
        for i in Seqs:
            print(i)

# --------------------------------------------------------------------


def calc_nbases(DNA):
    dna = Seq(DNA)
    if 'n' in dna or 'N' in dna:
        nbases = dna.count('n')+dna.count('N')
        print(nbases)

# --------------------------------------------------------------------


def filter_nbases(seq):
    seq = seq.replace("N", "")
    return seq
# --------------------------------------------------------------------


def reverse_complement(dnasequence):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N': 'N', 'n': 'n'}
    listofdnaletters = list(dnasequence)
    listofdnaletters = [basecomplement[base] for base in listofdnaletters]
    listofdnaletters = listofdnaletters[::-1]
    return ''.join(listofdnaletters)
# --------------------------------------------------------------------


def convert_To_fasta(filePath):
    if os.path.splitext(filePath)[1] == ".gbk":
        with open(filePath) as input_handle, open("{}.fasta".format(os.path.splitext(filePath)[0]), "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "genbank")
            count = SeqIO.write(sequences, output_handle, "fasta")
        print("converted %i records" % count)
    else:
        print("wrong path")

# --------------------------------------------------------------------


def online_alignment(seq, *FileName):
    try:
        fasta = seq
        result_handle = NCBIWWW.qblast('blastn', 'nt', fasta)
        blast_record = NCBIXML.read(result_handle)
        f = open(FileName[0], 'w')
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                f.write('\n')
                f.write('sequence: ')
                f.write(alignment.title)
                f.write('\n')
                f.write('length: ')
                f.write(str(alignment.length))
                f.write('\n')
                f.write('e value: ')
                f.write(str(hsp.expect))
                f.write('\n')
                f.write(hsp.query)
                f.write('\n')
                f.write(hsp.match)
                f.write('\n')
                f.write(hsp.sbjct)
        f.close()
    except:
        fasta = seq
        result_handle = NCBIWWW.qblast('blastn', 'nt', fasta)
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print("Sequence: ", alignment.title)
                print("/n")
                print("length: ", str(alignment.length))
                print("/n")
                print("e-value: ", str(hsp.expect))
                print("/n")
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)
# --------------------------------------------------------------------


argv = sys.argv[1:]
opts, args = getopt.gnu_getopt(argv, "o:")
if args[0] == 'transcribe':
    try:
        transcribe(args[1])
    except:
        print("transcribe function accept sequence as a parameter")
elif args[0] == 'gc' or args[0] == 'GC':
    try:
        print(gc(args[1]))
    except:
        print("GC function accept sequence as parameter")
elif args[0] == 'seq_alignment':
    try:
        try:
            if (opts[0][0] == '-o'):
                seq_alignment(args[1], args[2], opts[0][1])
                print("Alignment Done In The File")
        except:
            seq_alignment(args[1], args[2])
    except:
        print("seq_alignment function accept 2 sequences as parameters and you can add additional paramter by keyword -o")
elif args[0] == 'calc_nbases':
    try:
        calc_nbases(args[1])
    except:
        print("calc_nbases function accept sequence as a parameter")
elif args[0] == 'filter_nbases':
    try:
        print(filter_nbases(args[1]))
    except:
        print("filter_nbases function accept sequence as a parameter")
elif args[0] == 'convert_to_fasta':
    try:
        convert_To_fasta(args[1])
    except:
        print("convert_to_fasta function accept file as one parameter")
elif args[0] == 'is_valid':
    try:
        print(is_valid(args[1], args[2]))
    except:
        print("is_valid function accept sequence and type as a parameter")
elif args[0] == 'online_alignment':
    try:
        try:
            if (opts[0][0] == '-o'):
                online_alignment(args[1], opts[0][1])
                print("online_alignment DONE")
        except:
            online_alignment(args[1])
    except:
        print("online_alignment function take sequence parameter and optional paramater")
elif args[0] == 'reverse_complement':
    try:
        print(reverse_complement(args[1]))
    except:
        print("reverse_complement function accept sequence as a parameter")
elif args[0] == 'seq_alignment_files':
    try:
        try:
            if (opts[0][0] == '-o'):
                seq_alignment_files(args[1], args[2], opts[0][1])
                print("seq_alignment_files DONE")
        except:
            seq_alignment_files(args[1], args[2])
    except:
        print("seq_alignment_files accept 2 fasta files as a paramter and 1 output file as additional parameter")
elif args[0] == 'merge_fasta':
    try:
        try:
            if (opts[0][0] == '-o'):
                merge_fasta(args[1], args[2], opts[0][1], *args[3:])
                print("merge DONE")
        except:
            merge_fasta(args[1], args[2], "", *args[3:])
    except:
        print("merge_fasta functions accept 2 fasta files minumum as a paramter and 1 output file as additional parameter")
else:
    print("Function Not Found")
# print(args)
# print(opts)
# print(opts[0][1])
