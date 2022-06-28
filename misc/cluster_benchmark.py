import csv
from sklearn import metrics
import subprocess as sub
import argparse

parser = argparse.ArgumentParser(description='RATTLE running clustering and cluster summary step')
parser.add_argument('input', type=str, help='input reads file(required)')
parser.add_argument('output', type=str, help='output folder (required)')
parser.add_argument('threads', type=int, help='threads number to run RATTLE (required)')
parser.add_argument('--rna', action='store_true', help='whether to use RNA mode (instead of cDNA)')
args = parser.parse_args()

if args.rna:
    rattle_run = "./rattle cluster -i " + args.input + " -t " + str(args.threads) + " -o " + args.output + " --rna --iso"
else:
    rattle_run = "./rattle cluster -i " + args.input + " -t " + str(args.threads) + " -o " + args.output + " --iso"
p = sub.Popen(rattle_run, shell=True)
p.wait()
print("RATTLE isoform-level clustering completed")

rattle_run = "./rattle cluster_summary -i " + args.input + " -c " + args.output + "/clusters.out > " + args.output + "/summary.tsv"
p = sub.Popen(rattle_run, shell=True)
p.wait()
print("RATTLE cluster summary completed")

filename = args.output + "/summary.tsv"
csv_reader = csv.reader(open(filename))
tsp = {}
with open('./toyset/cluster_benchmark/input/ref.fa', 'r') as f:
    for count, line in enumerate(f, start=1):
        if count % 2 == 1:
            l = line.split()
            tsp[l[0][1:]] = l[3][5:]

labels_true_t = []
labels_pred =[]
labels_true_g = []
for line in csv_reader:
    labels_true_t.append(line[1])
    labels_pred.append(line[2])
    labels_true_g.append(tsp[line[1]])
print("homogeneity score with transcriptome is: {:.2f}%".format(metrics.homogeneity_score(labels_true_t, labels_pred) * 100))
print("completeness score with transcriptome is: {:.2f}%".format(metrics.completeness_score(labels_true_t, labels_pred) * 100))
print("homogeneity score with gene is: {:.2f}%".format(metrics.homogeneity_score(labels_true_g, labels_pred) * 100))
print("completeness score with gene is: {:.2f}%".format(metrics.completeness_score(labels_true_g, labels_pred) * 100))