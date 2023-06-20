def setup_samples(files, suffix):
    print('Processing sample information...')
    
    # samples setup

    df = pd.DataFrame(
        dict(id=files))
    df['batch'] = df['id'].str.replace('-', '_').str.split('_').apply(lambda x: x[0])
    df['sample_id'] = df['id'].apply(lambda x: x[: x.rfind('_')])
    df['dir'] = df['batch'].apply(illumina_path.get)
    df['path'] = df.apply(lambda x: f"{x.dir}/{x.id}_{suffix}", axis=1).iloc[0]

    df = df.set_index("sample_id").join(meta.set_index("sample_ID"))
    df['groups'] = df.loc[:, 'condition'].str.strip() + "_" + df.loc[:, 'replicate'].astype(int).astype(str)

    df.to_csv(
        "config/samples.tab",
        index=False
    )
