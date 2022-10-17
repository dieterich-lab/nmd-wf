

from os import symlink

df = pd.read_table(
    "config/samples.tab",
    sep='\t',
    index=False
)

illumina_id = df['id']
illumina_path = df['path']


rule symlink:
    input:
        bam=expand("{path}", path=illumina_id),
        bai=expand("{path}.bai", sample=sample, samples_dir=sample_path)
    output:
        bam=expand('mappings/{name}.bam', name=illumina_id),
        bai=expand('mappings/{name}.bam.bai', name=illumina_id)
    run:
        for bam_in, bai_in, bam_out, bai_out in zip(
                input.bam, input.bai, output.bam, output.bai):
            symlink(bam_in, bam_out)
            symlink(bai_in, bai_out)