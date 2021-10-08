include: "reference_alignment.smk"


rule gene_conversion_target_regions:
    input:
        genome=get_fai,
        bed=config["bed"],
    output:
        bed=temp("temp/{ref}/gene-conversion/target-regions.bed"),
    threads: 1
    conda:
        "../envs/env.yml"
    params:
        slop=2 * config.get("window", 10000),
    shell:
        """
        bedtools sort -i {input.bed} \
            | bedtools merge -i - \
            | bedtools slop -i - -g {input.genome} -b {params.slop} \
            | bedtools merge -i - \
            > {output.bed}
        """


rule make_query_windows:
    input:
        genome=get_fai,
        paf=rules.sam_to_paf.output.paf,
        bed=rules.gene_conversion_target_regions.output.bed,
    output:
        paf=temp("temp/{ref}/gene-conversion/{sm}_liftover.paf"),
    threads: 1
    conda:
        "../envs/env.yml"
    params:
        window=config.get("window", 10000),
        slide=config.get("slide", 5000),
    shell:
        """
        rb liftover --bed {input.bed} {input.paf} \
            | csvtk cut  -tT -f 1,3,4 \
            | bedtools makewindows -s {params.slide} -w {params.window} -b - \
            | rb liftover -q --bed /dev/stdin --largest {input.paf} \
            | grep -v "cg:Z:10000=" \
            > {output.paf}
        """


rule unzip_ref:
    input:
        query=get_asm,
    output:
        temp("temp/{ref}/gene-conversion/ref/{sm}_ref.fasta"),
    benchmark:
        "logs/{ref}/gene-conversion/alignment.{ref}_{sm}.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        seqtk seq -A -l 60 {input.query} > {output}
        samtools faidx {output}
        """


rule window_alignment:
    input:
        #ref=get_ref,
        ref=rules.alignment_index.output.mmi,
        query=rules.unzip_ref.output,
        paf=rules.make_query_windows.output.paf,
    output:
        aln=temp("temp/{ref}/gene-conversion/{sm}_windows.paf"),
    benchmark:
        "logs/{ref}/gene-conversion/alignment.{ref}_{sm}.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        minimap2 -K 1G -t {threads} \
            -cx asm20 \
            --secondary=no --eqx \
            {input.ref} \
                <( bedtools getfasta -name+ \
                    -fi {input.query} \
                    -bed <(awk -v OFS=$'\t' '{{name=$1":"$3"-"$4}}{{print $6,$8,$9,name}}' {input.paf}) \
                ) \
            > {output.aln}
        """


rule window_stats:
    input:
        paf=rules.window_alignment.output.aln,
        liftover_paf=rules.make_query_windows.output.paf,
    output:
        tbl="results/{ref}/gene-conversion/{sm}_windows.tbl.gz",
        liftover_tbl="results/{ref}/gene-conversion/{sm}_liftover_windows.tbl.gz",
    conda:
        "../envs/env.yml"
    threads: 4
    shell:
        """
        rb stats --paf {input.paf} \
            | pigz -p 4  > {output.tbl}
        rb stats --paf {input.liftover_paf} \
            | pigz -p {threads} > {output.liftover_tbl}
        """


rule candidate_gene_conversion:
    input:
        window=rules.window_stats.output.tbl,
        liftover=rules.window_stats.output.liftover_tbl,
    output:
        tbl=temp("temp/{ref}/gene-conversion/{sm}_candidate_windows.tbl"),
    conda:
        "../envs/env.yml"
    script:
        "../scripts/combine-mappings.R"


rule large_table:
    input:
        tbls=expand(
            rules.candidate_gene_conversion.output, sm=df.index, allow_missing=True
        ),
    output:
        tbl="results/{ref}/gene-conversion/all_candidate_windows.tbl.gz",
    conda:
        "../envs/env.yml"
    threads: 4
    shell:
        """
        (head -n 1 {input.tbls[0]}; tail -q -n +2 {input.tbls}) \
            | pigz -p {threads} \
            > {output.tbl}
        """


rule gene_conversion_windows:
    input:
        tbl=rules.large_table.output.tbl,
    output:
        tbl="results/{ref}/gene-conversion/gene_conversion_windows.tbl",
        interact="results/{ref}/gene-conversion/gene_conversion_interactions.bed",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/gene-conversion-windows.R"


rule make_big_bed:
    input:
        tbl=rules.gene_conversion_windows.output.tbl,
        interact=rules.gene_conversion_windows.output.interact,
        fai=get_fai,
    output:
        interact="results/{ref}/gene-conversion/all_candidate_interactions.bb",
        bb="results/{ref}/gene-conversion/all_candidate_windows.bb",
        bg="results/{ref}/gene-conversion/all_candidate_windows.bg",
        bw="results/{ref}/gene-conversion/all_candidate_windows.bw",
        bed=temp("temp/{ref}/gene-conversion/all_candidate_windows.bed"),
    conda:
        "../envs/env.yml"
    params:
        fmt=workflow.source_path("../scripts/bed.as"),
        interact=workflow.source_path("../scripts/interact.as"),
    threads: 4
    shell:
        """
        # make interactions
        bedtools sort -i {input.interact} \
            | awk '$3-$2 < 30e6' \
            > {output.bed}
        bedToBigBed -as={params.interact} \
            -type=bed5+13 {output.bed} {input.fai} {output.interact}

        # make others
        grep -v "reference_name" {input.tbl} \
            | bedtools sort -i - > {output.bed}

        bedToBigBed -as={params.fmt} -type=bed3+10 \
            {output.bed} {input.fai} {output.bb} 

        bedtools genomecov -i {output.bed} -g {input.fai} -bg > {output.bg}
        bedGraphToBigWig {output.bg} {input.fai} {output.bw}
        """


rule gene_conversion:
    input:
        expand(rules.large_table.output.tbl, ref=config.get("ref").keys()),
        expand(rules.gene_conversion_windows.output.tbl, ref=config.get("ref").keys()),
        expand(rules.make_big_bed.output, ref=config.get("ref").keys()),
    message:
        "Gene conversion run complete"
