
rule bed_to_pdf:
    input:
        bed=os.path.join(OUTPUT_DIR, "{ref}/bed/{sm}.bed"),
        script=workflow.source_path("scripts", "ideogram.R"),
        chm13_ktype=workflow.source_path("scripts", "chm13.karyo.RData"),
    output:
        pdf=os.path.join(OUTPUT_DIR, "{ref}/pdf/ideogram.{sm}.pdf"),
    threads: 1
    conda:
        "../envs/r.yml"
    shell:
        """
        Rscript {input.script} \
          --asm {input.bed} \
          --karyotype {input.chm13_ktype} \
          --plot {output.pdf}
        """
