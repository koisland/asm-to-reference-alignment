

rule alignment:
    input:
        ref=get_ref,
        query=get_asm,
    output:
        aln="temp/{ref}/{sm}.bam",
    log:
        "logs/alignment.{ref}_{sm}.log",
    benchmark:
        "benchmarks/alignment.{ref}_{sm}.benchmark.tsv"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    params:
        mm2_opts=config.get("mm2_opts", "-x asm20 --secondary=no -s 25000 -K 8G"),
    shell:
        """
        {{ minimap2 -t {threads} -a --eqx --cs \
            {params.mm2_opts} \
            {input.ref} {input.query} \
            | samtools view -F 4 -b -;}} \
            > {output.aln} 2> {log}
        """


rule alignment2:
    input:
        ref_fasta=get_ref,
        query=get_asm,
        aln=rules.alignment.output.aln,
    output:
        aln="temp/{ref}/{sm}.2.bam",
    log:
        "logs/alignment.{ref}_{sm}.2.log",
    benchmark:
        "benchmarks/alignment.{ref}_{sm}.2.benchmark.tsv"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    params:
        mm2_opts=config.get("mm2_opts", "-x asm20 --secondary=no -s 25000 -K 8G"),
        second_aln=config.get("second_aln", "no"),
    shell:
        """
        if [ {params.second_aln} == "yes" ]; then
          {{ minimap2 -t {threads} -a --eqx --cs \
              {params.mm2_opts} \
              <(seqtk seq \
                  -M <(samtools view -h {input.aln} | paftools.js sam2paf - | cut -f 6,8,9 | bedtools sort -i -) \
                  -n "N" {input.ref_fasta} \
              ) \
              <(seqtk seq \
                  -M <(samtools view -h {input.aln} | paftools.js sam2paf - | cut -f 1,3,4 | bedtools sort -i -) \
                  -n "N" {input.query} \
              ) \
              | samtools view -F 4 -b -;}}\
              > {output.aln} 2> {log}
        else
          samtools view -b -H {input.aln} > {output.aln}
        fi
        """


rule compress_sam:
    input:
        aln=rules.alignment.output.aln,
        aln2=rules.alignment2.output.aln,
    output:
        aln="results/{ref}/bam/{sm}.bam",
    threads: 1  # dont increase this, it will break things randomly 
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools cat {input.aln} {input.aln2} \
                 -o {output.aln}
        """
        # for some reason if I sort some cigars are turned into M instead of =/X
        #| samtools sort -m 8G --write-index \


rule sam_to_paf:
    input:
        aln=rules.compress_sam.output.aln,
    output:
        paf="results/{ref}/paf/{sm}.paf",
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools view -h {input.aln} \
            | paftools.js sam2paf - \
        > {output.paf}
        """


rule trim_and_break_paf:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        paf="results/{ref}/paf_trim_and_break/{sm}.paf",
    conda:
        "../envs/env.yml"
    params:
        break_paf=config.get("break_paf", 10_000),
    shell:
        """
        rustybam trim-paf {input.paf} \
            | rustybam break-paf --max-size {params.break_paf} \
        > {output.paf}
        """


rule aln_to_bed:
    input:
        aln=rules.compress_sam.output.aln,
    output:
        bed="results/{ref}/bed/{sm}.bed",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        rb --threads {threads} stats {input.aln} > {output.bed}
        """
