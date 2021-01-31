import re

import pandas as pd

RE_GENE_NAME = re.compile(r'Name=(?P<gene_name>.+?);')


def extract_gene_name(attributes_str):
    res = RE_GENE_NAME.search(attributes_str)
    return res.group('gene_name')


RE_GENE_ID = re.compile(r'gene_id=(?P<gene_id>ENSG.+?);')


def extract_gene_id(attributes_str):
    res = RE_GENE_ID.search(attributes_str)
    return res.group('gene_id')


RE_TRANS_GENE_ID = re.compile(r'Parent=gene:(?P<gene_id>ENSG.+?);')


def extract_transcript_gene_id(attributes_str):
    res = RE_TRANS_GENE_ID.search(attributes_str)
    return res.group('gene_id')


RE_EXON_TRANS_ID = re.compile(r'Parent=transcript:(?P<transcript_id>ENST\d+)')


def extract_exon_transcript_id(attributes_str):
    res = RE_EXON_TRANS_ID.search(attributes_str)
    return res.group('transcript_id')


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
    # print(df.head())
    # informacoes do dataframe
    # print(df.info())
    # primeira coluna seqid, dados da coluna de um dataframe
    df_unique = df.seqid.unique()
    # print(df_unique)
    # sequideos unicos, KI e GL sao sequencias de DNA - andaimes - no genoma que nao foram montadas com sucesso no
    # genoma
    # print(df_unique.shape)
    # contagem de variaveis categoricas
    # print(df.source.value_counts())
    # info sobre cromossomo inteiro
    gdf = df[df.source == 'GRCh38']
    # cada gdf corresponde a um dos 194 sequideos exclusivos
    # print(gdf.shape)
    # seleciona 10 entradas aleatorias
    # print(gdf.sample(10))
    gdf = gdf.copy()
    gdf['length'] = gdf.end - gdf.start + 1
    # 5 primeiras entradas
    # print(gdf.head())
    # print(gdf.length.sum())
    chrs = [str(_) for _ in range(1, 23)] + ['X', 'Y', 'MT']
    frsequnmounted = gdf[-gdf.seqid.isin(chrs)].length.sum() / gdf.length.sum()
    # comprimento total de entradas por entrada
    # fracao de sequencias nao montadas
    # print('Incompleto => ', frsequnmounted)
    edf = df[df.source.isin(['ensembl', 'havana', 'ensembl_havana'])]
    # print('edf sample: ', edf.sample(10))
    # elementos subgenicos
    # print('GENE:\n ', edf.type.value_counts())
    ndf = edf[edf.type == 'gene']

    ndf = ndf.copy()
    # print(ndf.sample(10).attributes.values)
    ndf['gene_name'] = ndf.attributes.apply(extract_gene_name)
    ndf['gene_id'] = ndf.attributes.apply(extract_gene_id)
    ndf['desc'] = ndf.attributes.apply(extract_description)
    # # remove coluna attributes
    ndf.drop('attributes', axis=1, inplace=True)
    # print(ndf.head())
    # print(ndf.shape)
    # print(ndf.gene_id.unique().shape)
    # print(ndf.gene_name.unique().shape)
    # count_df = ndf.groupby('gene_name').count().iloc[:, 0].sort_values().iloc[::-1]
    # print(type(gene_id_count))
    # print(count_df[count_df > 1].shape)
    # print(count_df.shape)
    #
    # print(count_df[count_df > 1].shape[0] / count_df.shape[0])
    # print(ndf[ndf.gene_name == 'SCARNA20'])
    # # duracao gene tipico
    # ndf['length'] = ndf.end - ndf.start + 1
    # print(ndf.length.describe())
    # ndf.length.plot(kind='hist', bins=50, logy=True)
    # # plt.show()
    # print(ndf[ndf.length > 2e6].sort_values('length').iloc[::-1])
    # print(ndf.sort_values('length').head())
    # # distribuicao de genes nos cromossomos
    # ndf = ndf[ndf.seqid.isin(chrs)]
    # chr_gene_counts = ndf.groupby('seqid').count().iloc[:, 0].sort_values().iloc[::-1]
    # print(chr_gene_counts)
    # print(df[(df.type == 'gene') & (df.seqid == 'MT')])
    # gdf = gdf[gdf.seqid.isin(chrs)]
    # gdf.drop(['start', 'end', 'score', 'strand', 'phase', 'attributes'], axis=1, inplace=True)
    # print(gdf.sort_values('length').iloc[::-1])
    # cdf = chr_gene_counts.to_frame(name='gene_count').reset_index()
    # print(cdf.head(2))
    # merged = gdf.merge(cdf, on='seqid')
    # print(merged)
    # print(merged[['length', 'gene_count']].corr())

    print(f'\n{"-" * 80}\nTRASNCRIPT')
    trans = edf[edf.type == 'transcript']
    trans = trans.copy()
    print(trans.head())
    trans['gene_id'] = trans.attributes.apply(extract_transcript_gene_id)
    trans.drop('attributes', axis=1, inplace=True)
    print('HEAD:\n', trans.head())
    gene_id_count = trans.groupby('gene_id').count().iloc[:, 0].sort_values().iloc[::-1]
    print('\nCOUNT TRANSCRIPT HEAD 10:\n', gene_id_count.head(10))
    print('N GENES COM TRANSCRIÇÕES: ', trans.gene_id.unique().shape[0])
    print('N TRANSCRIÇÕES: ', trans.gene_id.shape[0])
    print('N GENES TRANSCRIÇÕES = 1: ', gene_id_count[gene_id_count == 1].shape[0])
    print('N GENES TRANSCRIÇÕES > 1: ', gene_id_count[gene_id_count > 1].shape[0])
    print('TRANSCRIÇÕES MIN: ', min(gene_id_count.values))
    print('TRANSCRIÇÕES MAX: ', max(gene_id_count.values))
    print('TRANSCRIÇÕES MEDIA: ', trans.gene_id.shape[0] / trans.gene_id.unique().shape[0])
    print('% DE TRANSCRIÇÕES > 1: ', gene_id_count[gene_id_count > 1].shape[0] / trans.gene_id.unique().shape[0])

    print(f'\n{"-" * 80}\nEXON')
    exon = edf[edf.type == 'exon']
    exon = exon.copy()
    exon['transcript_id'] = exon.attributes.apply(extract_exon_transcript_id)
    exon['length'] = exon.end - exon.start + 1
    exon.drop('attributes', axis=1, inplace=True)
    print('HEAD:\n', exon.head())
    exon_transcript_id_count = exon.groupby('transcript_id').count().iloc[:, 0].sort_values().iloc[::-1]
    print('\nCOUNT EXON HEAD 10:\n', exon_transcript_id_count.head(10))
    print('N TRANSCRICOES COM EXONS: ', exon.transcript_id.unique().shape[0])
    print('N EXONS: ', exon.transcript_id.shape[0])
    print('N MAX EXON: ', max(exon_transcript_id_count.values))
    print('N MIN EXON: ', min(exon_transcript_id_count.values))
    print('N MEDIO EXON: ', exon.transcript_id.shape[0] / exon.transcript_id.unique().shape[0])
    print('TAMANHO MAX: ', max(exon.length))
    print('TAMANHO MIN', min(exon.length))
    print('TAMANHO MEDIO', exon.length.sum() / exon.shape[0])

    print(f'\n{"-" * 80}\nCDS')
    cds = edf[edf.type == 'CDS']
    cds = cds.copy()
    cds['transcript_id'] = cds.attributes.apply(extract_exon_transcript_id)
    cds['length'] = cds.end - cds.start + 1
    cds.drop('attributes', axis=1, inplace=True)
    print('HEAD:\n', cds.head())
    cds_transcript_id_count = cds.groupby('transcript_id').count().iloc[:, 0].sort_values().iloc[::-1]
    print('\nCOUNT CDS HEAD 10', cds_transcript_id_count.head(10))
    print('N TRANSCRICOES COM CDS: ', cds.transcript_id.unique().shape[0])
    print('N CDS: ', cds.transcript_id.shape[0])
    print('N MAX CDS: ', max(cds_transcript_id_count.values))
    print('N MIN CDS: ', min(cds_transcript_id_count.values))
    print('N MEDIO CDS: ', cds.transcript_id.shape[0] / cds.transcript_id.unique().shape[0])
    print('TAMANHO MAX: ', max(cds.length))
    print('TAMANHO MIN', min(cds.length))
    print('TAMANHO MEDIO', cds.length.sum() / cds.shape[0])

    print(f'\n{"-" * 80}\nUTR')
    utr = edf[edf.type.str.contains('UTR', regex=False)]
    utr = utr.copy()
    utr['transcript_id'] = utr.attributes.apply(extract_exon_transcript_id)
    utr['length'] = utr.end - utr.start + 1
    utr.drop('attributes', axis=1, inplace=True)
    print('HEAD:\n', utr.head())
    utr_transcript_id_count = utr.groupby('transcript_id').count().iloc[:, 0].sort_values().iloc[::-1]
    print('\nCOUNT UTR HEAD 10', utr_transcript_id_count.head(10))
    print('N TRANSCRICOES COM UTR: ', utr.transcript_id.unique().shape[0])
    print('N UTR: ', utr.transcript_id.shape[0])
    print('N MAX UTR: ', max(utr_transcript_id_count.values))
    print('N MIN UTR: ', min(utr_transcript_id_count.values))
    print('N MEDIO UTR: ', utr.transcript_id.shape[0] / utr.transcript_id.unique().shape[0])
    print('TAMANHO MAX: ', max(utr.length))
    print('TAMANHO MIN', min(utr.length))
    print('TAMANHO MEDIO', utr.length.sum() / utr.shape[0])
