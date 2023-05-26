def get_mix(wc):
    fname = paired[wc.name]
    return [
        f"mappings_selected/illumina/{fname[0]}.bam",
        f"mappings_selected/Nanopore/{fname[1]}.bam",
    ]

def get_illumina_stringtie(wc):
    files = Path("mappings/illumina/").glob("*.bam")
    files = [Path(x).resolve() for x in files]
    files = [str(x).replace("/Aligned.noS.bam", "_StringTie/stringtieB_outfile.gtf") for x in files]
    return files

def get_groups(wc):
    return ["mappings/illumina/{}_{}.bam".format(*x) 
        for x in illumina_mapping[wc.group]]

def extract_samples_replicates(samples, pattern=re.compile("^(.+)[_+$](.+)")):
    """
    Extract pairs of condition and replicate name from sample files
    :param str _pattern: pattern to . Default uses {condition}_{replicate} template
    :param list samples:
    :return:
    :rtype: list
    """
    return list(
        zip(*[re.match(
            pattern, x).groups() for x in samples]))


def get_raw_fastq(wc):
    r1 = (
        df.pipe(pd.DataFrame.set_index, 'CCG_Sample_ID')
          .pipe(lambda x: x.loc[wc.name, "raw_reads"]))
    r1 = str(r1)
    r2 = r1.replace("_R1_", "_R2_")
            
    return {"r1": r1, "r2": r2}




