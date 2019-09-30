import sys
import re
from tqdm import tqdm
import mmap
import argparse

parser = argparse.ArgumentParser(description='Calculate known/novel splice sites from PAF alignment and ref GTF file')

parser.add_argument('ref_gtf', type=str,
                    help='Reference GTF file (required)')

parser.add_argument('aln_paf', type=str,
                    help='Alignment file in PAF format (required)')

parser.add_argument('--beautiful', action='store_true',
                    help='Beautiful output (instead of csv)')

args = parser.parse_args()

txtExons = {} # txtExons[chr][tid] = [(x,y), (x2,y2)...]
totalExonsInRef = 0

knownExons = {}
knownIntrons = {}
knownTranscriptsIntronLevel = {}
knownTranscriptsExonLevel = {}

pafExons = {}
pafIntrons = {}
pafTranscriptsIntronLevel = {}
pafTranscriptsExonLevel = {}

singleExons = 0

if args.beautiful:
    print("")

def get_num_lines(file_path):
    with open(file_path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# read gtf reference
with open(args.ref_gtf) as gtf:
    pbar = tqdm(gtf, total=get_num_lines(args.ref_gtf))
    pbar.set_description_str("[INFO] Parsing GTF...")

    for line in pbar:
        if line[0] == "#":
            continue

        info = line.split()    

        if info[2] == "exon":
            tid = info[11].replace('"', "").replace(";", "")
            
            chrom = info[0]
            start = int(info[3]) - 1
            end = int(info[4]) - 1
            strand = info[5]

            if chrom not in txtExons:
                txtExons[chrom] = {}
            
            if tid not in txtExons[chrom]:
                txtExons[chrom][tid] = []
            
            txtExons[chrom][tid].append((start,end))
            totalExonsInRef += 1

with tqdm(total = totalExonsInRef) as pbar:
    pbar.set_description_str("[INFO] Processing exons...")

    for chrom in txtExons:
        for tid in txtExons[chrom]:
            tidExons = chrom
            tidIntrons = chrom
            lastExonEnd = -1

            txtExons[chrom][tid].sort(key=lambda x: x[0])
            for exon in txtExons[chrom][tid]:
                start = exon[0]
                end = exon[1]

                knownExons[chrom + "," + str(start) + "-" + str(end)] = True
                tidExons += "," + str(start) + "-" + str(end)

                # intron level
                if lastExonEnd != -1:
                    intronStart = lastExonEnd
                    intronEnd = start - 1
                    lastExonEnd = end + 1

                    knownIntrons[chrom + "," + str(intronStart) + "-" + str(intronEnd)] = True
                    tidIntrons += "," + str(intronStart) + "-" + str(intronEnd)
                else:
                    lastExonEnd = end + 1

            # add chains
            if tidExons in knownTranscriptsExonLevel:
                tqdm.write("[WARN] Detected two transcripts with the same exon chain: " + knownTranscriptsExonLevel[tidExons] + " -> " + tid, file=sys.stderr)
            knownTranscriptsExonLevel[tidExons] = tid

            if tidIntrons != chrom:
                if tidIntrons in knownTranscriptsIntronLevel:
                    tqdm.write("[INFO] Detected two transcripts with the same intron chain: " + knownTranscriptsIntronLevel[tidIntrons] + " -> " + tid, file=sys.stderr)
                knownTranscriptsIntronLevel[tidIntrons] = tid
            else:
                singleExons += 1

            pbar.update(len(txtExons[chrom][tid]))
                

tqdm.write("[INFO] Single-exon transcripts: " + str(singleExons), file=sys.stderr)

# read paf aln
with open(args.aln_paf) as paf:
    pbar = tqdm(paf, total=get_num_lines(args.aln_paf))
    pbar.set_description_str("[INFO] Parsing PAF...")

    for line in pbar:
        info = line.split()

        # build exons from cigar string
        chrom = info[5]
        strand = info[4]
        start = int(info[7])
        end = start - 1

        cigar = info[-1].split("cg:Z:")[1].split()[0]
        cigarInfo = [(int(v[0]), v[1]) for v in re.findall(r"([0-9]+)([A-Z=]+)", cigar, re.I)]
        readExons = {}
        readExonsChain = ""
        readIntrons = {}
        readIntronsChain = ""

        for pair in cigarInfo:
            if pair[1] == "M":
                end += pair[0]
            elif pair[1] == "D":
                end += pair[0]
            elif pair[1] == "N":
                k = chrom + "," + str(start) + "-" + str(end)
                readExons[k] = 1

                if readExonsChain == "":
                    readExonsChain = k
                else:
                    readExonsChain += "," + str(start) + "-" + str(end)

                intronStart = end + 1
                intronEnd = intronStart + pair[0] - 1
                ki = chrom + "," + str(intronStart) + "-" + str(intronEnd)
                readIntrons[ki] = 1

                if readIntronsChain == "":
                    readIntronsChain = ki
                else:
                    readIntronsChain += "," + str(intronStart) + "-" + str(intronEnd)

                start = end + pair[0] + 1
                end = start - 1
            else:
                if pair[1] != "I":
                    tqdm.write("[ERR] Unsupported CIGAR op " + str(pair[0]) + pair[1], file=sys.stderr)

        # final exon
        if end - start > 1:
            k = chrom + "," + str(start) + "-" + str(end)
            readExons[k] = 1

            if readExonsChain == "":
                readExonsChain = k
            else:
                readExonsChain += "," + str(start) + "-" + str(end)

        for k in readExons:
            if k in pafExons:
                pafExons[k] += 1
            else:
                pafExons[k] = 1
        
        if readExonsChain in pafTranscriptsExonLevel:
            pafTranscriptsExonLevel[readExonsChain] += 1
        else:
            pafTranscriptsExonLevel[readExonsChain] = 1

        for k in readIntrons:
            if k in pafIntrons:
                pafIntrons[k] += 1
            else:
                pafIntrons[k] = 1
        
        if readIntronsChain != "":
            if readIntronsChain in pafTranscriptsIntronLevel:
                pafTranscriptsIntronLevel[readIntronsChain] += 1
            else:
                pafTranscriptsIntronLevel[readIntronsChain] = 1

        # if info[0] == "5c5b4bb4-42e0-4c9f-bac4-a4c2e59c8b6a":
        #     print(readIntrons)
        
        # if info[0] == "8459ffb7-9211-4063-92be-150030a9fd7b/6cc575bdf77decb4e6fec8129022e38229f87aa4":
        #     print(readExons)

if not args.beautiful:
    print("level,known_in_ref,unique_in_reads,ref_found,p_ref_found,total_in_reads,known_in_total_reads,novel_in_total_reads,p_known_in_total_reads,p_novel_in_total_reads")

pafKnown = {}
pafUnknown = {}
pafCountKnown = 0
pafCountUnknown = 0
for e in pafIntrons:
    if e in knownIntrons:
        pafKnown[e] = True
        pafCountKnown += pafIntrons[e]
    else:
        pafUnknown[e] = True
        pafCountUnknown += pafIntrons[e]

if args.beautiful:
    print("")
    print("########################################")
    print("#             INTRON LEVEL             #")
    print("########################################")
    print("Introns in reference: " + str(len(knownIntrons)))
    print("Unique introns in reads: " + str(len(pafIntrons)))
    print("Reference introns found: {:d}/{:d} ({:.2f}%)".format(len(pafKnown), len(knownIntrons), float(len(pafKnown)) * 100.0 / float(len(knownIntrons))))
    # print("Novel exons found: {:d}/{:d} ({:.2f}%)".format(len(pafUnknown), len(pafExons), float(len(pafUnknown)) * 100.0 / float(len(pafExons))))
    print("Total introns in reads: " + str(pafCountKnown + pafCountUnknown))
    print("--> Known: {:d} ({:.2f}%)".format(pafCountKnown, float(pafCountKnown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
    print("--> Novel: {:d} ({:.2f}%)".format(pafCountUnknown, float(pafCountUnknown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
else:
    csvLine = "intron"
    csvLine += ",{:d}".format(len(knownIntrons))
    csvLine += ",{:d}".format(len(pafIntrons))
    csvLine += ",{:d}".format(len(pafKnown))
    csvLine += ",{:.2f}".format(float(len(pafKnown)) / float(len(knownIntrons)))
    csvLine += ",{:d}".format(pafCountKnown + pafCountUnknown)
    csvLine += ",{:d}".format(pafCountKnown)
    csvLine += ",{:d}".format(pafCountUnknown)
    csvLine += ",{:.2f}".format(float(pafCountKnown) / float(pafCountKnown + pafCountUnknown))
    csvLine += ",{:.2f}".format(float(pafCountUnknown) / float(pafCountKnown + pafCountUnknown))
    print(csvLine)

pafKnown = {}
pafUnknown = {}
pafCountKnown = 0
pafCountUnknown = 0
for e in pafExons:
    if e in knownExons:
        pafKnown[e] = True
        pafCountKnown += pafExons[e]
    else:
        pafUnknown[e] = True
        pafCountUnknown += pafExons[e]

if args.beautiful:
    print("\n")
    print("########################################")
    print("#              EXON LEVEL              #")
    print("########################################")
    print("Exons in reference: " + str(len(knownExons)))
    print("Unique exons in reads: " + str(len(pafExons)))
    print("Reference exons found: {:d}/{:d} ({:.2f}%)".format(len(pafKnown), len(knownExons), float(len(pafKnown)) * 100.0 / float(len(knownExons))))
    # print("Novel exons found: {:d}/{:d} ({:.2f}%)".format(len(pafUnknown), len(pafExons), float(len(pafUnknown)) * 100.0 / float(len(pafExons))))
    print("Total exons in reads: " + str(pafCountKnown + pafCountUnknown))
    print("--> Known: {:d} ({:.2f}%)".format(pafCountKnown, float(pafCountKnown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
    print("--> Novel: {:d} ({:.2f}%)".format(pafCountUnknown, float(pafCountUnknown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
else:
    csvLine = "exon"
    csvLine += ",{:d}".format(len(knownExons))
    csvLine += ",{:d}".format(len(pafExons))
    csvLine += ",{:d}".format(len(pafKnown))
    csvLine += ",{:.2f}".format(float(len(pafKnown)) / float(len(knownExons)))
    csvLine += ",{:d}".format(pafCountKnown + pafCountUnknown)
    csvLine += ",{:d}".format(pafCountKnown)
    csvLine += ",{:d}".format(pafCountUnknown)
    csvLine += ",{:.2f}".format(float(pafCountKnown) / float(pafCountKnown + pafCountUnknown))
    csvLine += ",{:.2f}".format(float(pafCountUnknown) / float(pafCountKnown + pafCountUnknown))
    print(csvLine)



pafKnown = {}
pafUnknown = {}
pafCountKnown = 0
pafCountUnknown = 0
for e in pafTranscriptsIntronLevel:
    if e in knownTranscriptsIntronLevel:
        pafKnown[e] = True
        pafCountKnown += pafTranscriptsIntronLevel[e]
    else:
        pafUnknown[e] = True
        pafCountUnknown += pafTranscriptsIntronLevel[e]

if args.beautiful:
    print("\n")
    print("########################################")
    print("#   TRANSCRIPT LEVEL (INTRON CHAIN)    #")
    print("########################################")
    print("Transcripts in reference: " + str(len(knownTranscriptsIntronLevel)))
    print("Unique transcripts in reads: " + str(len(pafTranscriptsIntronLevel)))
    print("Reference transcripts found: {:d}/{:d} ({:.2f}%)".format(len(pafKnown), len(knownTranscriptsIntronLevel), float(len(pafKnown)) * 100.0 / float(len(knownTranscriptsIntronLevel))))
    print("Total transcripts in reads: " + str(pafCountKnown + pafCountUnknown))
    print("--> Known: {:d} ({:.2f}%)".format(pafCountKnown, float(pafCountKnown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
    print("--> Novel: {:d} ({:.2f}%)".format(pafCountUnknown, float(pafCountUnknown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
else:
    csvLine = "intron_chain"
    csvLine += ",{:d}".format(len(knownTranscriptsIntronLevel))
    csvLine += ",{:d}".format(len(pafTranscriptsIntronLevel))
    csvLine += ",{:d}".format(len(pafKnown))
    csvLine += ",{:.2f}".format(float(len(pafKnown)) / float(len(knownTranscriptsIntronLevel)))
    csvLine += ",{:d}".format(pafCountKnown + pafCountUnknown)
    csvLine += ",{:d}".format(pafCountKnown)
    csvLine += ",{:d}".format(pafCountUnknown)
    csvLine += ",{:.2f}".format(float(pafCountKnown) / float(pafCountKnown + pafCountUnknown))
    csvLine += ",{:.2f}".format(float(pafCountUnknown) / float(pafCountKnown + pafCountUnknown))
    print(csvLine)


pafKnown = {}
pafUnknown = {}
pafCountKnown = 0
pafCountUnknown = 0
for e in pafTranscriptsExonLevel:
    if e in knownTranscriptsExonLevel:
        pafKnown[e] = True
        pafCountKnown += pafTranscriptsExonLevel[e]
    else:
        pafUnknown[e] = True
        pafCountUnknown += pafTranscriptsExonLevel[e]

if args.beautiful:
    print("\n")
    print("########################################")
    print("#    TRANSCRIPT LEVEL (EXON CHAIN)     #")
    print("########################################")
    print("Transcripts in reference: " + str(len(knownTranscriptsExonLevel)))
    print("Unique transcripts in reads: " + str(len(pafTranscriptsExonLevel)))
    print("Reference transcripts found: {:d}/{:d} ({:.2f}%)".format(len(pafKnown), len(knownTranscriptsExonLevel), float(len(pafKnown)) * 100.0 / float(len(knownTranscriptsExonLevel))))
    print("Total transcripts in reads: " + str(pafCountKnown + pafCountUnknown))
    print("--> Known: {:d} ({:.2f}%)".format(pafCountKnown, float(pafCountKnown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
    print("--> Novel: {:d} ({:.2f}%)".format(pafCountUnknown, float(pafCountUnknown) * 100.0 / float(pafCountKnown + pafCountUnknown)))
    print("")
else:
    csvLine = "exon_chain"
    csvLine += ",{:d}".format(len(knownTranscriptsExonLevel))
    csvLine += ",{:d}".format(len(pafTranscriptsExonLevel))
    csvLine += ",{:d}".format(len(pafKnown))
    csvLine += ",{:.2f}".format(float(len(pafKnown)) / float(len(knownTranscriptsExonLevel)))
    csvLine += ",{:d}".format(pafCountKnown + pafCountUnknown)
    csvLine += ",{:d}".format(pafCountKnown)
    csvLine += ",{:d}".format(pafCountUnknown)
    csvLine += ",{:.2f}".format(float(pafCountKnown) / float(pafCountKnown + pafCountUnknown))
    csvLine += ",{:.2f}".format(float(pafCountUnknown) / float(pafCountKnown + pafCountUnknown))
    print(csvLine)