import numpy as np
import pandas as pd
import os
import shutil
import docker
from Bio import motifs
import time

### search for motifs in promoter sequences of coregulated genes; search for discovered motifs in all promoter sequences
def find_motifs(initial_GRN, gene_promoters):
    print('Motifs search...')
    genes = initial_GRN.index
    n_of_reg_genes = np.sum(initial_GRN != 0)
    verified_GRN_np = np.zeros((np.shape(initial_GRN)[0], np.shape(initial_GRN)[1]), dtype = int)
    verified_GRN = pd.DataFrame(verified_GRN_np, columns=genes, index=genes)
    client = docker.from_env()
    client.images.pull('memesuite/memesuite:latest')
    sto_file = open('output/discovered_motifs.sto', 'w') # motif sequences in Stockholm format

    for r in range (0, np.shape(initial_GRN)[0]):
        if n_of_reg_genes[r] >= 5: # motif search only if TFs regulates at least 5 genes
            #print(f'Motif: {r+1} / {np.shape(initial_GRN)[0]}')
            ofile = open('temporary_coreg_seq.fasta', 'w') # promoter sequences
            prom_no = 1
            coreg_seqs = 0
            for c in range(0, np.shape(initial_GRN)[1]): # write individual promoters into temporary_coreg_seq.fasta, filter out promoters shorter than 5
                if initial_GRN.iloc[r, c] > 0 and len(gene_promoters[c]) >= 5:
                    ofile.write('>' + str(genes[r]) + '_' + str(prom_no) + '\n' + str(gene_promoters[c]) + '\n')
                    coreg_seqs += 1
                    prom_no += 1
            ofile.close()
            # motif discovery using MEME Suite Docker if promoters were exported to the fasta file
            if (os.path.getsize('temporary_coreg_seq.fasta') > 0) & (coreg_seqs > 1): # number of sequences must be 2 at least
                if not os.path.exists('meme_out'):
                    os.mkdir('meme_out')
                os.chmod('meme_out', 0o777)
                os.chmod('temporary_coreg_seq.fasta', 0o777)
                client.containers.run('memesuite/memesuite:latest',
                                                  'meme -mod zoops -minw 5  -dna temporary_coreg_seq.fasta',
                                                  volumes={os.getcwd(): {'bind': '/home/meme', 'mode': 'rw'}}, detach=True)
                sleep = 0
                while not os.path.exists('meme_out/meme.xml'):
                    time.sleep(2)
                    sleep += 1
                    if sleep > 30: # skip current motif search as meme suite has not returned results in a long time
                        break
                try:
                    handle = open('meme_out/meme.xml')
                    record = motifs.parse(handle, 'meme')
                    handle.close()
                except:
                    continue

                for m in range(0, len(record)):  # search each discovered motif
                    motif = record[m]
                    motif.pseudocounts = 0.5
                    pssm = motif.pssm # Position-Specific Scoring Matrix
                    distribution = pssm.distribution()
                    threshold = distribution.threshold_fpr(0.01)

                    for n in range(0, len(gene_promoters)): # search motif in each promoter
                        if len(gene_promoters[n]) >= motif.instances[0].length:
                            score = pssm.calculate(gene_promoters[n])
                            try:
                                score_max = max(score)
                            except TypeError: # score has only 1 value
                                score_max = score
                            except ValueError: # score has 0 values
                                score_max = threshold - 10
                            rpssm = pssm.reverse_complement()
                            score_rev = rpssm.calculate(gene_promoters[n])
                            try:
                                score_rev_max = max(score_rev)
                            except TypeError:
                                score_rev_max = score_rev
                            except ValueError:
                                score_rev_max = threshold - 10
                            if score_max >= threshold or score_rev_max >= threshold:
                                verified_GRN.iloc[r, n] = 1 # assign found motif in the promoter to the TF

                # MEME Suite output into Stockholm file format
                motifs_stockholm(genes, r, record, sto_file)
                try:
                    shutil.rmtree('meme_out')
                except PermissionError:
                    pass
            os.remove('temporary_coreg_seq.fasta')
    sto_file.close()
    print('Motifs search done.')
    return(verified_GRN)

### MEME Suite output into Stockholm file format
def motifs_stockholm(genes, r, record, sto_file):
    for i in range(0, len(record)):  # save each discovered motif into Stockholm file
        motifs_discovered = record[i]
        sto_file.write('# STOCKHOLM 1.0' + '\n')
        sto_file.write('#=GF ID   ' + str(genes[r]) + '\n')

        for z in range(0,motifs_discovered.num_occurrences): # save each sequence containing discovered motif into Stockholm file
            sto_file.write(
                str(motifs_discovered.instances[z].sequence_name) + '	' + str(motifs_discovered.instances[z]) + '\n')
        sto_file.write('//' + '\n')
