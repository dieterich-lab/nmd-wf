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

def bam_to_fastq(wc):

    # file = Path(str(df_illumina[df_illumina.group =l= wc.illumina_name]['files']))
    name = group2sample[wc.illumina_name]
    print(name)
    file = Path(f'/prj/Niels_Gehring/nmd_transcriptome/mappings/illumina/{name}.bam')
    file = file.resolve()
    file = file.with_name('Log.out')
    with open(file) as fin:
        for line in fin:
            if line.startswith("##### Command Line:"):
                line = next(fin)
                break
    for args in line.split("--"):
        if args.startswith("readFilesIn"):
            break
    args = args.split()        
    
    return {'r1': file.parent.parent.parent / args[1], 
            'r2': file.parent.parent.parent / args[2]}


def bam_to_raw_fastq(wc):
    r1, r2 = bam_to_fastq(wc)
