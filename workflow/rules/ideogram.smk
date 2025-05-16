
rule bed_to_pdf:
    input:
        bed=os.path.join(OUTPUT_DIR, "{ref}/bed/{sm}.bed"),
        script=workflow.source_path(os.path.join("..", "scripts", "ideogram.R")),
        chm13_ktype=workflow.source_path(os.path.join("..", "scripts", "chm13.karyo.RData")),
    output:
        pdf=os.path.join(OUTPUT_DIR, "{ref}/pdf/ideogram.{sm}.pdf"),
    params:
        ideogram_min=config.get("ideogram_min", 1e6)
    threads: 1
    conda:
        "../envs/r.yml"
    shell:
        """
        Rscript {input.script} \
          --asm {input.bed} \
          --karyotype {input.chm13_ktype} \
          --min {params.ideogram_min} \
          --plot {output.pdf}
        """
