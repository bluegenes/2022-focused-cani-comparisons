import os
import pandas as pd
# compare different ANI methods on a series of genomes

out_dir = "output.ani-compare"
logs_dir = os.path.join(out_dir, "logs")

comparison_file = "gtdb-comparisons.largest-size-diff.n10.csv"
comparisons = pd.read_csv(comparison_file)
identA= comparisons['identA'].tolist()
identB= comparisons['identB'].tolist()
IDENTS = list(set(identA+identB))

original_genomedir = "/home/ntpierce/2021-rank-compare/genbank/genomes"

# read in comparison file; get all IDENTS + use them to copy genoms + run pyani/fastani
# parse pyani, fastani output for desired comparisons
# report
# smash: use signatures from k21 zipfile? or just rebuild sigs from fasta files?

rule all:
    input: 
        os.path.join(out_dir, "genomes", "genome-filepaths.txt"),
        expand(os.path.join(out_dir, "genomes", "{acc}_genomic.fna"), acc=IDENTS),
        os.path.join(out_dir, "pyani", "ANIb_results/pyani.csv"),
        os.path.join(out_dir, "fastani","fastani.tsv"),
        os.path.join(out_dir, "sourmash", "signatures.zip"),

def get_all_genomes(w):
    comparison_accs = IDENTS
    genomes = []
    for acc in comparison_accs:
        genome_file = os.path.join(original_genomedir, f"{acc}_genomic.fna.gz")
        genomes.append(genome_file)
    return genomes

### split into genome folders + unzip fna files and generate classes/labels, then run pyANI index and compare
# make a folder with all genomes
rule copy_and_unzip_genomes:
    input:
        get_all_genomes
    output:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt"),
        genomes=expand(os.path.join(out_dir, "genomes", "{acc}_genomic.fna"), acc=IDENTS),
    params:
        acc_list = IDENTS,
        genome_dir = lambda w: os.path.join(out_dir, 'genomes'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2",
    run:
        import hashlib
        os.makedirs(params.genome_dir, exist_ok=True)
        with open(output.classes, 'w') as out_classes:
            with open(output.labels, 'w') as out_labels:
                for acc in params.acc_list:
                    fn = os.path.join(original_genomedir, f"{acc}_genomic.fna.gz") 
                    dest = os.path.join(params.genome_dir, f"{acc}_genomic.fna.gz")
                    dest_unz = os.path.join(params.genome_dir, f"{acc}_genomic.fna")
                    md5_file = os.path.join(params.genome_dir, f"{acc}_genomic.md5")
                    shell("cp {fn} {dest}")
                    shell("gunzip -c {fn} > {dest_unz}")
                    # get md5 of unzipped fna
                    with open(dest_unz, "rb") as f:
                        bytes = f.read()
                        md5 = hashlib.md5(bytes).hexdigest()
                    # write to md5 file
                    with open(md5_file, 'w') as mfile:
                        mfile.write(f"{md5}\t{dest_unz}\n")
                    fna_base = os.path.basename(dest_unz).rsplit('.fna')[0]
                    out_classes.write(f"{md5}\t{fna_base}\t{acc}\n")
                    out_labels.write(f"{md5}\t{fna_base}\t{acc}\n")

# we copy all the fasta files to the genome_dir for pyani, so we can just use those paths here
rule build_genome_filepaths:
    output: os.path.join(out_dir, "genomes", "genome-filepaths.txt")
    params:
        genome_dir=os.path.join(out_dir, "genomes")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2",
    run:
        with open(str(output), "w") as out:
            acc_list = IDENTS
            for acc in acc_list:
                fn = os.path.join(params.genome_dir, f"{acc}_genome.fna")
                out.write(f"{fn}\n")

# sourmash
# rebuild sigs so we can use k21, k31, k51, scaled = 1 so we can do scaled= 1,10,1000 comparisons

rule sourmash_sketch:
    input: os.path.join(out_dir, "genomes", "genome-filepaths.txt")
    output: os.path.join(out_dir, "sourmash", "signatures.zip")
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, "sourmash", "sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash", "sketch.benchmark")
    shell:
        """
        sourmash sketch fromfile {input} -p k=21,k=31,k=51,scaled=1 -o {output} > {log}
        """

# use api so we can get all values
rule sourmash_api_compare:
    input: 
        sigs=os.path.join(out_dir, "sourmash", "signatures.zip"),
        comparisons=comparison_file,#COMPARISON_FILES
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, "sourmash", "api-compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "api-compare.benchmark")
    shell:
        """
        python sourmash-api-compare.py {input.sigs} --comparisons {input.comparisons} -o {output} 2> {log}
        """

# just use python api -- not very many comparisons
#rule sourmash_compare_using_api:
#    input:

rule pyani_index_and_createdb:
    input:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt")
    output:
        classes=os.path.join(out_dir, "genomes", "classes.txt"),
        labels=os.path.join(out_dir, "genomes", "labels.txt"),
        db=os.path.join(out_dir, "genomes", ".pyanidb")
    params:
        genome_dir = lambda w: os.path.join(out_dir, 'genomes'),
        pyanidb = lambda w: os.path.join(out_dir, 'genomes',f".pyanidb"),
        classes_basename = "classes.txt",
        labels_basename = "labels.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2",
    log: os.path.join(logs_dir, "pyani", "index-and-createdb.log")
    benchmark: os.path.join(logs_dir, "pyani", "index-and-createdb.benchmark")
    conda: "conf/env/pyani0.3.yml"
    shell:
        """
        pyani index -i {params.genome_dir} --classes {params.classes_basename} --labels {params.labels_basename}
        pyani createdb --dbpath {params.pyanidb} -v -l {log}
        """

rule pyANI_ANIb:
    input:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt"),
        idx_classes=os.path.join(out_dir, "genomes", "classes.txt"),
        idx_labels=os.path.join(out_dir, "genomes", "labels.txt"),
    output:
        covF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_similarity_errors.tab"),
        bn =  os.path.join(out_dir, "pyani/ANIb_results","blastn_output.tar.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        time=1200,
        partition="med2"
    params:
        #pyanidb = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, f".pyani-{w.path}/pyanidb"),
        genome_dir = lambda w: os.path.join(out_dir, 'genomes'),
        output_dir = lambda w: os.path.join(out_dir, 'pyani', "ANIb_results"),
    log: os.path.join(logs_dir, "pyani_anib", "pyANI-anib.log")
    benchmark: os.path.join(logs_dir, "pyani_anib", "pyANI-anib.benchmark")
    conda: "conf/env/pyani0.2.yml"
    shell:
        """
        average_nucleotide_identity.py -i {params.genome_dir} \
             -o {params.output_dir} -m ANIb -v \
             --labels {input.labels} --classes {input.classes} \
             --force > {log}
        """

#localrules: aggregate_ranktax_anib
rule aggregate_ranktax_anib:
    input:
        covF= os.path.join(out_dir, "pyani/ANIb_results", "ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_similarity_errors.tab"),
    output:
        os.path.join(out_dir, "pyani", "ANIb_results/pyani.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/ANIb_results"),
        #compare_rank = compare_rank,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2"
    shell:
        """
        python aggregate-pyani-results.py {params.results_dir} 
        """
#--rank {params.compare_rank} --ranktax-name {wildcards.ranktax} --output-csv {output} --pyani-version v0.2

rule compare_via_fastANI:
    input: os.path.join(out_dir, "genomes", "genome-filepaths.txt")
    output: os.path.join(out_dir, "fastani","fastani.tsv"),
    threads: 11
    resources:
        mem_mb=lambda wildcards, attempt: attempt *9000,
        time=12000,
        partition="bmm",
    log: os.path.join(logs_dir, "fastani", "fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "fastani.benchmark")
    conda: "conf/env/fastani.yml"
    shell:
        """
        fastANI --threads {threads} --ql {input} --rl {input} -o {output} 2> {log}
        """


#localrules: parse_fastani_ranktax
#rule parse_fastani_ranktax:
#    input: os.path.join(out_dir, "fastani", "fastani.tsv")
#    output: os.path.join(out_dir, "fastani", "fastani.parsed.csv")
#    #params:
#        #compare_rank= compare_rank,
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *3000,
#        time=12000,
#        partition="med2",
#    log: os.path.join(logs_dir, "fastani", "parse_fastani.log")
#    benchmark: os.path.join(logs_dir, "fastani", "parse_fastani.benchmark")
#    shell:
#        """
#        python aggregate-fastani.py --fastani-ranktax {input} --rank {params.compare_rank} \
#                                    --ranktax-name {wildcards.ranktax} \
#                                    --output-csv {output} 2> {log}
#        """

