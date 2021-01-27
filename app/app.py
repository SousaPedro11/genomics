import re

import matplotlib.pyplot as plt
import pandas as pd

RE_GENE_NAME = re.compile(r'Name=(?P<gene_name>.+?);')


def extract_gene_name(attributes_str):
    res = RE_GENE_NAME.search(attributes_str)
    return res.group('gene_name')


RE_GENE_ID = re.compile(r'gene_id=(?P<gene_id>ENSG.+?);')


def extract_gene_id(attributes_str):
    res = RE_GENE_ID.search(attributes_str)
    return res.group('gene_id')


RE_DESC = re.compile('description=(?P<desc>.+?);')


def extract_description(attributes_str):
    res = RE_DESC.search(attributes_str)
    if res is None:
        return ''
    else:
        return res.group('desc')


def read_file():
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
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
    return df


def run_genomics():
    df = read_file()
    # 5 primeiras entradas
    print(df.head())
    # informacoes do dataframe
    print(df.info())
    # primeira coluna seqid, dados da coluna de um dataframe
    df_unique = df.seqid.unique()
    print(df_unique)
    # sequideos unicos, KI e GL sao sequencias de DNA - andaimes - no genoma que nao foram montadas com sucesso no
    # genoma
    print(df_unique.shape)
    # contagem de variaveis categoricas
    print(df.source.value_counts())
    # info sobre cromossomo inteiro
    gdf = df[df.source == 'GRCh38']
    # cada gdf corresponde a um dos 194 sequideos exclusivos
    print(gdf.shape)
    # seleciona 10 entradas aleatorias
    print(gdf.sample(10))
    gdf = gdf.copy()
    gdf['length'] = gdf.end - gdf.start + 1
    # 5 primeiras entradas
    print(gdf.head())
    print(gdf.length.sum())
    chrs = [str(_) for _ in range(1, 23)] + ['X', 'Y', 'MT']
    frsequnmounted = gdf[-gdf.seqid.isin(chrs)].length.sum() / gdf.length.sum()
    # comprimento total de entradas por entrada
    # fracao de sequencias nao montadas
    print('Incompleto => ', frsequnmounted)
    edf = df[df.source.isin(['ensembl', 'havana', 'ensembl_havana'])]
    print(edf.sample(10))
    # elementos subgenicos
    print('GENE => ', edf.type.value_counts())
    ndf = edf[edf.type == 'gene']
    ndf = ndf.copy()
    print(ndf.sample(10).attributes.values)
    ndf['gene_name'] = ndf.attributes.apply(extract_gene_name)
    ndf['gene_id'] = ndf.attributes.apply(extract_gene_id)
    ndf['desc'] = ndf.attributes.apply(extract_description)
    # remove coluna attributes
    ndf.drop('attributes', axis=1, inplace=True)
    print(ndf.head())
    print(ndf.shape)
    print(ndf.gene_id.unique().shape)
    print(ndf.gene_name.unique().shape)
    count_df = ndf.groupby('gene_name').count().iloc[:, 0].sort_values().iloc[::-1]
    print(count_df.head(10))
    print(count_df[count_df > 1].shape)
    print(count_df.shape)
    print(count_df[count_df > 1].shape[0] / count_df.shape[0])
    print(ndf[ndf.gene_name == 'SCARNA20'])
    # duracao gene tipico
    ndf['length'] = ndf.end - ndf.start + 1
    print(ndf.length.describe())
    ndf.length.plot(kind='hist', bins=50, logy=True)
    # plt.show()
    print(ndf[ndf.length > 2e6].sort_values('length').iloc[::-1])
    print(ndf.sort_values('length').head())
    # distribuicao de genes nos cromossomos
    ndf = ndf[ndf.seqid.isin(chrs)]
    chr_gene_counts = ndf.groupby('seqid').count().iloc[:, 0].sort_values().iloc[::-1]
    print(chr_gene_counts)
    print(df[(df.type == 'gene') & (df.seqid == 'MT')])
    gdf = gdf[gdf.seqid.isin(chrs)]
    gdf.drop(['start', 'end', 'score', 'strand', 'phase', 'attributes'], axis=1, inplace=True)
    print(gdf.sort_values('length').iloc[::-1])
    cdf = chr_gene_counts.to_frame(name='gene_count').reset_index()
    print(cdf.head(2))
    merged = gdf.merge(cdf, on='seqid')
    print(merged)
    print(merged[['length', 'gene_count']].corr())
