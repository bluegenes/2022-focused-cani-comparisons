import os
import pandas as pd
# compare different ANI methods on a series of genomes

out_dir = "output.ani-compare-test"
logs_dir = os.path.join(out_dir, "logs")

#comparison_file = "gtdb-comparisons.largest-size-diff.n10.csv"
comparison_file = "family.test-comparisons.n5.csv"
comparisons = pd.read_csv(comparison_file)
basename= "comparisons"
identA= comparisons['identA'].tolist()
identB= comparisons['identB'].tolist()
#IDENTS = identA+identB
IDENTS = list(dict.fromkeys(identA + identB)) # rm duplicates but preserve order to make sure copy_genomes rule doesn't get rerun
KSIZE = [21]
SCALED = [1,1000]


original_genomedir = "/home/ntpierce/2021-rank-compare/genbank/genomes"

# read in comparison file; get all IDENTS + use them to copy genoms + run pyani/fastani
# parse pyani, fastani output for desired comparisons
# report
# smash: use signatures from k21 zipfile? or just rebuild sigs from fasta files?

rule all:
    input: 
        #os.path.join(out_dir, "genomes", "genome-filepaths.txt"),
        #expand(os.path.join(out_dir, "genomes", "{acc}_genomic.fna"), acc=IDENTS),
        #os.path.join(out_dir, "sourmash", "signatures.zip"),
        expand(os.path.join(out_dir, "sourmash", f"{basename}.k{{ksize}}-sc{{scaled}}.cANI.csv"), ksize=[21], scaled=[1,1000]),
        os.path.join(out_dir, "combinedANI.csv"),

### split into genome folders + unzip fna files and generate classes/labels, then run pyANI index and compare
# make a folder with all genomes
rule copy_and_unzip_genomes:
    input:
        expand(os.path.join(original_genomedir, "{acc}_genomic.fna.gz"),acc=IDENTS)
    output:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt"),
        genomes=expand(os.path.join(out_dir, "genomes", "{acc}_genomic.fna"), acc=IDENTS),
    params:
        acc_list = IDENTS,
        genome_dir = os.path.join(out_dir, 'genomes'),
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
                    shell("cp {fn} {dest}") # do we need both gzipped and unzipped versions??
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
    output:
        filepaths=os.path.join(out_dir, "genomes", "genome-filepaths.txt"),
        fromfile_csv=os.path.join(out_dir, "genomes", "fromfile.csv")
    params:
        genome_dir=os.path.join(out_dir, "genomes")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2",
    run:
        with open(str(output.filepaths), "w") as out, open(str(output.fromfile_csv), 'w') as ff_out:
            acc_list = IDENTS
            ff_out.write('name,genome_filename,protein_filename\n')
            for acc in acc_list:
                fn = os.path.join(params.genome_dir, f"{acc}_genomic.fna")
                out.write(f"{fn}\n")
                ff_out.write(f"{acc},{fn},\n")

# sourmash
# rebuild sigs so we can use k21, k31, k51, scaled = 1 so we can do scaled= 1,10,1000 comparisons

rule sourmash_sketch:
    input: os.path.join(out_dir, "genomes", "fromfile.csv")
    output: os.path.join(out_dir, "sourmash", "signatures.zip")
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, "sourmash", "sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash", "sketch.benchmark")
    shell:
        """
        sourmash sketch fromfile {input} -p dna,k=21,k=31,k=51,scaled=1 -o {output} > {log}
        """

# just use python api -- not very many comparisons
rule sourmash_api_compare:
    input: 
        sigs=os.path.join(out_dir, "sourmash", "signatures.zip"),
        comparisons=comparison_file,#COMPARISON_FILES
    output: 
        compare_csv=os.path.join(out_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.comparison-info.csv"),
        ani_csv=os.path.join(out_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.cANI.csv")
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.api-compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.api-compare.benchmark")
    shell:
        """
        python sourmash-api-compare.py {input.sigs} -c {input.comparisons} -o {output.compare_csv} --ani-csv {output.ani_csv} 2> {log}
        """

rule pyani_index_and_createdb:
    input:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt")
    output:
        classes=os.path.join(out_dir, "genomes", "classes.txt"),
        labels=os.path.join(out_dir, "genomes", "labels.txt"),
        db=os.path.join(out_dir, "genomes", ".pyanidb"),
    params:
        genome_dir = os.path.join(out_dir, 'genomes'),
        pyanidb = os.path.join(out_dir, 'genomes',f".pyanidb"),
        classes_basename = "classes.txt",
        labels_basename = "labels.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        time=1200,
        partition="bmm",
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
        comparisons=comparison_file,
        covF= os.path.join(out_dir, "pyani/ANIb_results", "ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_similarity_errors.tab"),
    output:
        os.path.join(out_dir, "pyani", "pyani.ANIb.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/ANIb_results"),
        #compare_rank = compare_rank,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2"
    shell:
        """
        python aggregate-pyani-results.py {params.results_dir} --pyani-version 0.2 -c {input.comparisons} -o {output}
        """
#--rank {params.compare_rank} --ranktax-name {wildcards.ranktax} --output-csv {output} --pyani-version v0.2

rule compare_via_fastANI:
    input: os.path.join(out_dir, "genomes", "genome-filepaths.txt")
    output: os.path.join(out_dir, "fastani","fastani.tsv"),
    threads: 11
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        time=12000,
        partition="bmm",
    log: os.path.join(logs_dir, "fastani", "fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "fastani.benchmark")
    conda: "conf/env/fastani.yml"
    shell:
        """
        fastANI --threads {threads} --ql {input} --rl {input} -o {output} 2> {log}
        """


localrules: parse_fastani
rule parse_fastani:
    input: 
        fastani=os.path.join(out_dir, "fastani", "fastani.tsv"),
        comparisons=comparison_file,
    output: os.path.join(out_dir, "fastani", "fastani.ANI.csv")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=12000,
        partition="med2",
    log: os.path.join(logs_dir, "fastani", "parse_fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "parse_fastani.benchmark")
    shell:
        """
        python parse-fastani-results.py --fastani {input.fastani} -c {input.comparisons} \
                                        --output-csv {output} 2> {log}
        """

localrules: combine_ani
rule combine_ani:
    input: 
        fastani= os.path.join(out_dir, "fastani", "fastani.ANI.csv"),
        pyani=os.path.join(out_dir, "pyani", "pyani.ANIb.csv"),
        sourmash=expand(os.path.join(out_dir, "sourmash", f"{basename}.k{{ksize}}-sc{{scaled}}.cANI.csv"), ksize=KSIZE, scaled=SCALED),
    output: os.path.join(out_dir, "combinedANI.csv")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=12000,
        partition="med2",
    log: os.path.join(logs_dir, "combine-ani", "combine-ani.log")
    benchmark: os.path.join(logs_dir, "combine-ani", "combine-ani.benchmark")
    shell:
        """
        python combine-ani.py --fastani {input.fastani} --pyani {input.pyani} \
                              --sourmash {input.sourmash} --output-csv {output} 2> {log}
        """

