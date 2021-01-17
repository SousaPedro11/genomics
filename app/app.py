import pandas as pd


def read_file():
    print(pd.__version__)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    col_names = [
        'seqid',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase',
        'attributes'
    ]

    df = pd.read_csv(
        'Homo_sapiens.GRCh38.85.gff3.gz',
        compression='gzip',
        sep='\t',
        comment='#',
        low_memory=False,
        header=None,
        names=col_names
    )
    print(df.head())
