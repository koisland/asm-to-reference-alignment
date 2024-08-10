

rule alignment:
    input:
        ref=lambda wc: REFERENCES[str(wc.ref)],
        query=lambda wc: SAMPLES[str(wc.sm)],
    output:
        aln=os.path.join(TEMP_DIR, "{ref}", "{sm}.bam"),
    log:
        os.path.join(LOGS_DIR, "alignment_{ref}_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARK_DIR, "alignment_{ref}_{sm}.tsv")
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


rule bam_to_paf:
    input:
        aln=rules.alignment.output.aln,
    output:
        paf=os.path.join(OUTPUT_DIR, "{ref}/paf/{sm}.paf"),
    log:
        os.path.join(LOGS_DIR, "bam_to_paf_{ref}_{sm}.log"),
    conda:
        "../envs/env.yml"
    shell:
        """
        {{ samtools view -h {input.aln} | paftools.js sam2paf - ;}} > {output.paf} 2> {log}
        """


rule trim_and_break_paf:
    input:
        paf=rules.bam_to_paf.output.paf,
    output:
        paf=os.path.join(OUTPUT_DIR, "{ref}/paf_trim_and_break/{sm}.paf"),
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOGS_DIR, "trim_and_break_paf_{ref}_{sm}.log"),
    params:
        break_paf=config.get("break_paf", 10_000),
    shell:
        """
        {{ rb trim-paf {input.paf} \
            | rb break-paf --max-size {params.break_paf} \
        ;}} > {output.paf} 2> {log}
        """


rule aln_to_bed:
    input:
        aln=rules.alignment.output.aln,
    output:
        bed=os.path.join(OUTPUT_DIR, "{ref}/bed/{sm}.bed"),
    log:
        os.path.join(LOGS_DIR, "aln_to_bed_{ref}_{sm}.log"),
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        rb --threads {threads} stats {input.aln} > {output.bed} 2> {log}
        """
