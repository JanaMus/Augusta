from numpy import sum, zeros, shape, max
from pandas import DataFrame
from os import path, mkdir, chmod, getcwd, remove
from shutil import rmtree
from docker import from_env
from Bio import motifs
from time import sleep
from gc import collect

### search for motifs in promoter sequences of coregulated genes; search for discovered motifs in all promoter sequences
def find_motifs(initial_GRN, gene_promoters, motifs_max_time):
    print('Motifs search...')
    genes = initial_GRN.index
    n_of_reg_genes = sum(initial_GRN != 0)
    verified_GRN_np = zeros((shape(initial_GRN)[0], shape(initial_GRN)[1]), dtype='int8')
    verified_GRN = DataFrame(verified_GRN_np, columns=genes, index=genes, dtype='int8')
    client = from_env()
    client.images.pull('memesuite/memesuite:latest')
    sto_file = open('output/discovered_motifs.sto', 'w') # motif sequences in Stockholm format

    for r in range (0, shape(initial_GRN)[0]):
        if n_of_reg_genes[r] >= 5: # motif search only if TFs regulates at least 5 genes
            gene_name = genes[r]
            #print(f'Motif: {r+1} / {shape(initial_GRN)[0]}')
            with open('temporary_coreg_seq.fasta', 'w') as ofile: # promoter sequences
                prom_no = 1
                coreg_seqs = 0
                for c in range(0, shape(initial_GRN)[1]): # write individual promoters into temporary_coreg_seq.fasta, filter out promoters shorter than 5
                    if initial_GRN.iloc[r, c] > 0 and len(gene_promoters[c]) >= 5:
                        ofile.write('>' + str(gene_name) + '_' + str(prom_no) + '\n' + str(gene_promoters[c]) + '\n')
                        coreg_seqs += 1
                        prom_no += 1

            # motif discovery using MEME Suite Docker if promoters were exported to the fasta file
            if (path.getsize('temporary_coreg_seq.fasta') > 0) & (coreg_seqs > 1): # number of sequences must be 2 at least
                if not path.exists('meme_out'):
                    mkdir('meme_out')
                chmod('meme_out', 0o777)
                chmod('temporary_coreg_seq.fasta', 0o777)
                container = client.containers.run('memesuite/memesuite:latest',
                                                  'meme -mod zoops -minw 5  -dna temporary_coreg_seq.fasta',
                                                  volumes={getcwd(): {'bind': '/home/meme', 'mode': 'rw'}},
                                                  detach=True, name=str(gene_name), remove=True)

                sleep_count = 0
                while client.containers.list(filters={'name':str(gene_name)}): # check if MEME Suite is still running
                    if sleep_count > motifs_max_time:  # skip current motif search as meme suite container has not finished in a long time
                        print("Motif search terminated for TF: ", str(gene_name), " (time limit exceeded).")
                        try:
                            container.stop(timeout=5)
                        except: # container already finished or removal already in progress
                            pass
                        break
                    else:
                        sleep(5)
                        sleep_count += 5
                del container

                try:
                    with open('meme_out/meme.xml') as handle:
                        record = motifs.parse(handle, 'meme')
                except:
                    continue

                for m in range(0, len(record)):  # search each discovered motif
                    motif = record[m]
                    motif.pseudocounts = 0.5
                    pssm = motif.pssm # Position-Specific Scoring Matrix
                    distribution = pssm.distribution()
                    threshold = round(distribution.threshold_balanced(1000),2) # distribution.threshold_fpr(0.01) - puvodne, pak balanced(10000) ###smazat KOMENTAR
                    threshold_int = int(threshold*100)

                    for n in range(0, len(gene_promoters)): # search motif in each promoter
                        if len(gene_promoters[n]) >= motif.instances[0].length:
                            score = pssm.calculate(gene_promoters[n])
                            score_max = round(max(score, initial=threshold),2)
                            score_max_int = int(score_max*100)
                            rpssm = pssm.reverse_complement() # reverse complement
                            score_rev = rpssm.calculate(gene_promoters[n])
                            score_rev_max = round(max(score_rev, initial=threshold),2)
                            score_rev_max_int = int(score_rev_max*100)
                            if (score_max_int > threshold_int) or (score_rev_max_int > threshold_int):
                                verified_GRN.iloc[r, n] = 1 # assign found motif in the promoter to the TF

                    # MEME Suite output into Stockholm file format
                    motifs_stockholm(gene_name, m, motif.consensus, sto_file)
                try:
                    rmtree('meme_out')
                except PermissionError:
                    pass
        collect()
    remove('temporary_coreg_seq.fasta')
    sto_file.close()
    print('Motifs search done.')
    return(verified_GRN)

### MEME Suite output into Stockholm file format
def motifs_stockholm(gene_name, count, consensus, sto_file): # save each discovered motif into Stockholm file
    sto_file.write('# STOCKHOLM 1.0' + '\n')                # header
    sto_file.write('#=GF ID   ' + str(gene_name) + '\n')    # header
    sto_file.write(                                         # sequence name, consensus motif
                str(gene_name) + '_' + str(count) +
                '	' + str(consensus) + '\n')
    sto_file.write('//' + '\n')                             # tail
