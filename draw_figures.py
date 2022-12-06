"""
1. Heatmap of NKT subtype markers.
2. Swearm plots of Klf2, Zeb2 by clusters
3. Correlation comparison with Zeb2 or Klrg1
4. Scatter plot of scRNA-seq coodinates
"""

import os, sys, re
import scipy.io
import scipy.sparse
import json
import pathlib
import plotly
import plotly.subplots
import plotly.graph_objs as go
import plotly.express as px
import argparse
import pandas as pd
import numpy as np

import scipy.stats
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

sys.path.append('/mnt/nas/genomeplatform/scripts')

import tkutil
import cellrangerwrapper
    

def swarmplot_by_cluster():#filename_expr, clusters, **kwargs):
    """
params: 
    filename_expr: count or normalized value for expression
    clusters : cluster name and column indexes
    auto_open : flag to open graph
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--seurat-dir', default='merged/seurat', metavar='directory', help='seurat.R script output directory')
    parser.add_argument('-c', default=None, metavar='cluster file', help='precalculated cluster fiile in tsv format')
    parser.add_argument('--cluster', nargs='+', default=['C7', 'C11', 'C1', 'C9', 'C10', 'C2'])
    parser.add_argument('-e', default=None, metavar='tsv file', help='Gene expression table (TSV format)')
    parser.add_argument('--verbose', action='store_true', help='verbosicty')
    parser.add_argument('--open', action='store_true', help='Open HTML after drawing a chart')    
    parser.add_argument('-o', default='out', metavar='directory', help='output directory')
    parser.add_argument('-x', type=int, default=20, metavar='number', help='number of graphs')
    parser.add_argument('-g', nargs='+', default=['Klf2'], metavar='string', help='Genes to show graphs')
    parser.add_argument('--column', default=2, type=int, metavar='number', help='number of columns of subplots')
    parser.add_argument('--mode', choices=['strip', 'swarm'], default='strip', help='Graph mode (swarm/swtrip plot)')
    args, cmd = parser.parse_known_args()

    max_genes = args.x
    n_cols = args.column
    auto_open = args.open
    outdir = pathlib.Path(args.o)
    outdir.mkdir(exist_ok=True)
    seurat = pathlib.Path(args.seurat_dir)
    mode = args.mode

    genes = args.g
    if os.path.isfile(genes[0]):
        df = pd.read_csv(genes[0], sep='\t')
        for c_ in ('Gene', 'gene', 'genes', 'gene_id'):
            if c_ in df.columns:
                genes = list(df[c_].values.reshape(-1)[0:max_genes])
                break
    fn_clu = pathlib.Path(args.c or (seurat / 'seurat.clusters.tsv').as_posix())
    fn_expr = pathlib.Path(args.e or (seurat / 'seurat.normalized.tsv').as_posix())
    fn_log = outdir / f'run.{mode}.log'
    fn_chart = outdir / f'{mode}.html'
    fn_pdf = outdir / f'{mode}.pdf'
    fn_info = outdir / f'run.{mode}.info'
    logger = tkutil.get_logger(os.path.basename(__file__), fn_log)
    if args.verbose: logger.setLevel(10)
    logger.info(', '.join(genes))

    info = {
        'command':sys.argv,
        'genes':genes,
        'data': {
            'cluster':fn_clu.as_posix(),
            'expression':fn_expr.as_posix(),
        },
        'chart':{
            'columns':n_cols,
            'filename':fn_chart.as_posix(),
            'pdf':fn_pdf.as_posix(),
        }
    }

    clu = pd.read_csv(fn_clu, sep='\t')
    cluval = clu.values.reshape(-1)
    clu2idx = {}
    for i in range(clu.shape[0]):
        c = cluval[i]
        if c not in clu2idx: clu2idx[c] = []
        clu2idx[c].append(i)
    logger.info('{} clusters of {} cells detected'.format(len(clu2idx), clu.shape[0]))

    accepted_clusters = []
    n_cells = 0
    n_clusters = 0
    if args.cluster is None:
        clusters = [int(_) for _ in sorted(clu2idx.keys())]
    else:
        clusters = args.cluster

    for c in clusters:
        if re.match('C\\d+$', c):
            c = int(c[1:]) - 1
            accepted_clusters.append(c)
        elif c.isdigit():
            accepted.clusters.append(int(c))
        else:
            raise Exception('cluster format should be either digits or CX (X > 0)')
    expr = cellrangerwrapper.load_gene_expression(fn_expr, genes=genes)
    logger.info('{}x{} expression table loaded'.format(expr.shape[0], expr.shape[1]))
    n_clusters = len(clu2idx)
    

    titles = [g for g in genes if g in expr.index]
    n_rows = (len(titles) + n_cols - 1) // n_cols
    figs = plotly.subplots.make_subplots(cols=n_cols, rows=n_rows, subplot_titles=titles)

    col = row = 1
    pxtable = []
    for gene in titles:
        for i, c in enumerate(accepted_clusters):
            idx = clu2idx[c]
            values = expr.loc[gene].values.reshape(-1)
            values = values[idx]
            # print(values)
            n_cells += len(idx)
            # df = pd.DataFrame(np.array([[i] * len(idx), values]), index=['X', 'Y']).T
            # print(df)
            # f_ = px.strip(df, x='X', y='Y')
            # print(f_.keys())
            # figs.add_trace(f_)
            for v in values:
                pxtable.append([gene, 'C{}'.format(c+1), v])
            # figs.add_trace(px.strip(df, x='X', y='Y'))
            # figs.add_trace(px.box(x=[i] * len(idx), y=values, points='all'))
            # figs.add_trace(px.strip(x=[i] * len(idx), y=values, hover_name='C{}'.format(c + 1)), col=col, row=row)
            col += 1
            if col > n_cols:
                col = 1
                row += 1
    pxdf = pd.DataFrame(pxtable, columns=['Gene', 'Cluster', 'Count'])
    f = px.strip(data_frame=pxdf, x='Cluster', y='Count', log_y=True, labels='Gene', facet_col='Gene', facet_col_wrap=n_rows)
    f.update_layout(plot_bgcolor='white')
    f.update_xaxes(linewidth=1, linecolor='black', mirror='ticks')
    f.update_yaxes(linewidth=1, linecolor='black', mirror='ticks')
    plotly.offline.plot(f, filename=fn_chart.as_posix(), auto_open=auto_open)


    
    colors = 'navy', 'limegreen', 'gold', 'orange', 'brown', 'cyan'
    colors = ['navy', ]#'cadetblue', 'navy', ]
    linecolors = ['orange', ]

    # sns.set()
    # print(axes)
    # print(n_rows, n_cols)
    # print(axes.shape)

    # fn_pdf = outdir / f'swarm.pdf'
    sns.axes_style('white')
    sns.set(rc={'figure.figsize':(n_cols * 4 + 1, n_rows * 4+ 1), 'axes.facecolor':'white', 'axes.edgecolor':'black', 'axes.facecolor':'white'})
    doc = matplotlib.backends.backend_pdf.PdfPages(fn_pdf.as_posix())

    fig, axes = plt.subplots(n_rows, n_cols)
    if len(axes.shape) == 1:
        ax = [axes, ]
    else:
        ax = axes
    # ax_index = 0
    row = col = 0
    fig.subplots_adjust(hspace=1.0, wspace=0.5)
    for gene in titles:
        pxtable = []
        means = []

        for i, c in enumerate(accepted_clusters):
            idx = clu2idx[c]
            values = expr.loc[gene].values.reshape(-1)
            values = values[idx]
            # print(values)
            n_cells += len(idx)
            # df = pd.DataFrame(np.array([[i] * len(idx), values]), index=['X', 'Y']).T
            # print(df)
            # f_ = px.strip(df, x='X', y='Y')
            # print(f_.keys())
            # figs.add_trace(f_)
            for v in values:
                pxtable.append(['C{}'.format(c+1), v])
            means.append(np.mean(values))
            # figs.add_trace(px.strip(df, x='X', y='Y'))
            # figs.add_trace(px.box(x=[i] * len(idx), y=values, points='all'))
            # figs.add_trace(px.strip(x=[i] * len(idx), y=values, hover_name='C{}'.format(c + 1)), col=col, row=row)

        df = pd.DataFrame(pxtable, columns=['Cluster', 'Count'])
        vals_ = df['Count'].values.reshape(-1)
        pt_ = max(1, np.percentile(vals_, 99))
        yfig = int(np.log10(pt_) )
        ymax = 10 ** int(yfig)
        while ymax < pt_:
            for d in 1, 2, 5:
                ymax = d * (10 ** int(yfig))
                # print(ymax, pt_)
                if ymax > pt_:
                    break
            yfig += 1
        if mode == 'swarm':
            sampled = df.sample(1000).sort_index()
        # print(df.shape, sampled.shape)
            ax_ = sns.swarmplot(data=sampled, x='Cluster', y='Count', palette=colors, size=2,
            edgecolor='black', marker='o', ax=ax[row][col], zorder=1).set_title(gene)
        else:
            ax_ = sns.stripplot(data=df, x='Cluster', y='Count', palette=colors, size=1.5, jitter=0.25, dodge=True,
            edgecolor='black', marker='o', ax=ax[row][col], zorder=1).set_title(gene)
        # print(ax_)
        for i, c in enumerate(accepted_clusters):
            ax[row][col].hlines(means[i], i-.35, i+.35,  zorder=100, colors=linecolors)
        # df_mean = df.groupby('Cluster', sort=False)['Count'].mean()
        # sx_.hlines(y, i-.25, i+.25, zorder=2) for i, y in df_mean.reset_index()

        ax[row][col].set(ylim=[0,int(ymax)])
        # ax[row][col].annotate(gene, xy=(0,1.02), xytext=(x1,1.02), xycoords='axes fraction')



        ax[row][col].set_ylabel('Normalized count')
        ax[row][col].set_xlabel('Cluster')


        # ax_index += 1
        col += 1
        if col >= n_cols:
            col = 0
            row += 1

        # fn_pdf = outdir / f'swarm_{gene}.pdf'
        # df = pd.DataFrame(pxtable, columns=['Gene', 'Cluster', 'Count'])
        # with sns.axes_style('whitegrid'):
        #     ax = sns.stripplot(data=pxdf, x='Cluster', y='Count', palette=colors,
        #     edgecolor='black', marker='o', )
            # ax = sns.stripplot(data=pxdf, x='Cluster', y='Count', size=3.0, jitter=0.25, dodge=False, palette=colors,
            # edgecolor='black', marker='o', )
            # ax.set_yticks([0,.25, .5, .75, 1.0])
            # ax.set_ylabel('Normalized count')
            # ax.set_xlabel('Cluster')
            # print(ax.get_xticklabels())
            # ax.set_xticklabels(colidx)
            # for p0, p1 in ((0, 1), (2,3)):
            #     # vi = df[df['Sample']==colidx[p0]]['Distance'].values.reshape(-1)
            #     # vj = df[df['Sample']==colidx[p1]]['Distance'].values.reshape(-1)
            #     # pval = scipy.stats.brunnermunzel(vi, vj).pvalue
            #     # plt.plot([p0, p1], [1.02, 1.02], 'k-', lw=2)
            #     # ax.text(p0 + .25, 1.05, '{:.2e}'.format(pval))#1.5, 1.1, 'text1')
            #     x0 = (p0 + .5) / 4
            #     x1 = (p1 + .5) / 4
            #     ax.annotate('', xy=(x0,1.02), xytext=(x1,1.02), xycoords='axes fraction',
            #     arrowprops=dict(arrowstyle='-', color='black'))
    doc.savefig(pad_inches=1, bbox_inches='tight')
    doc.close()
    plt.clf()

    info['data']['cells'] = int(n_cells)
    info['data']['clusters'] = clusters
    # plotly.offline.plot(figs, filename=fn_out, auto_open=auto_open)

    # print(json.dumps(info, indent=2))

    with open(fn_info, 'w') as fo:
        json.dump(info, fo, indent=2)

def draw_correlation_map():
    """Draw scatter plot of given genes.
    parameters
    -e [filename] : expression tsv file
    -c [filename] : cluster tsv file
    -o [filename] : Chart HTML file
    -g gene1 gene2... : Genes
    --cluster cls1 cls2 ... : Clusters used for calculation
    """
    parser = argparse.ArgumentParser()        
    parser.add_argument('-e', default='merged/seurat/seurat.normalized.tsv')
    parser.add_argument('-c', default='merged/seurat/seurat.clusters.tsv')
    parser.add_argument('-o', default='correlation.html')
    parser.add_argument('-g', nargs='+', default=['Klf2', 'Ly6e', 'Zeb2', 'S1pr5', 'Ctla2a'])
    parser.add_argument('--cluster', nargs='+', default=['C7', 'C11', 'C1', 'C9', 'C10', 'C2'])
    args, cmd = parser.parse_known_args()

    expr = cellrangerwrapper.load_gene_expression(args.e, genes=args.g)
    available = []
    for g in args.g:
        if g not in expr.index or g in available: continue
        available.append(g)
    genes = available

    clu = pd.read_csv(args.c, sep='\t', index_col=0)
    genes = args.g
    n_genes = len(genes)
    n_charts = n_genes * (n_genes - 1) // 2
    n_rows = int(np.ceil(np.sqrt(n_charts)))
    n_cols = (n_charts + n_rows - 1) // n_rows

    cluval = clu.values.reshape(-1)
    clu2idx = {}
    accepted_clusters = []
    for i in range(clu.shape[0]):
        c = cluval[i]
        if c not in clu2idx: clu2idx[c] = []
        clu2idx[c].append(i)
    if args.cluster is None:
        clusters = [int(_) for _ in sorted(clu2idx.keys())]
    else:
        clusters = args.cluster
    for c in clusters:
        if re.match('C\\d+$', c):
            c = int(c[1:]) - 1
            accepted_clusters.append(c)
        elif c.isdigit():
            accepted.clusters.append(int(c))
        else:
            raise Exception('cluster format should be either digits or CX (X > 0)')

    figs = plotly.subplots.make_subplots(cols=n_cols, rows=n_rows)

    col = 1
    row = 1
    color = []
    idxlist = []
    colors = ['orange', 'cyan', 'lime', 'darkslateblue', 'brown', 'darkgreen', 'lightgray', 
        'magenta', 'cyan']
    for k,c in enumerate(accepted_clusters):
        color += [colors[k]] * len(clu2idx[c])
        idxlist.append(np.array(clu2idx[c], dtype=np.int32))
    idx = np.concatenate(idxlist)

    for i, g in enumerate(genes):
        vali = expr.loc[g].values[idx].reshape(-1)
        for j in range(i + 1, len(genes)):
            valj = expr.loc[genes[j]].values[idx].reshape(-1)
            figs.add_trace(go.Scatter(x=vali, y=valj, mode='markers', name='',
                marker=dict(size=3, color=color)), row=row, col=col)
            figs.update_xaxes(row=row, col=col, title=g, linewidth=1, linecolor='black', mirror='ticks')
            figs.update_yaxes(row=row, col=col, title=genes[j], linewidth=1, linecolor='black', mirror='ticks')
            r = scipy.stats.pearsonr(vali, valj)[0]
            # print(i, j, row, col, g, genes[j], r)
            col += 1
            if col > n_cols:
                col = 1
                row += 1
    figs.update_layout(plot_bgcolor='white', showlegend=False)
    plotly.offline.plot(figs, filename=args.o)

def _load_mt_counts(filename):
    import hashlib, gzip, pickle
    md5 = hashlib.md5()
    md5.update((os.path.abspath(filename) + ';count_and_mtcounts').encode('utf-8'))
    code = md5.hexdigest()[0:6]    
    fn_cache = os.path.join(os.path.dirname(filename), f'.{code}.cache')
    print(fn_cache)
    if os.path.exists(fn_cache) and os.path.getsize(fn_cache) > 1000:
        with gzip.open(fn_cache) as fi:
            df = pickle.load(fi)
        if df.shape[1] > 1000 and df.shape[0] >= 2:
            return df

    mtcounts = None
    nuccounts = None
    barcodes = None
    with open(filename) as fi:
        barcodes = fi.readline().strip().split('\t')
        n_cells = len(barcodes)
        mtcounts = [0] * n_cells
        nuccounts = [0] * n_cells
        for line in fi:
            items = line.split('\t')
            sys.stderr.write('\033[K{}\r'.format(items[0]))
            if items[0].strip('"').lower().startswith('mt-'):
                vals = mtcounts
            else:
                vals = nuccounts
            for i, val in enumerate(items[1:]):
                if val != '"0"' and i < n_cells:
                    vals[i] += int(val.strip('"'))
        sys.stderr.write('\033[K\r')
    df = pd.DataFrame([nuccounts, mtcounts], index=['mt', 'nuc'], columns=barcodes)
    with gzip.open(fn_cache, 'wb') as fo:
        pickle.dump(df, fo)
    return df


def draw_cell_state():
    """
    MT-ratio vs read count map to identify dead / apoptotic /duplex cells.
    -c [filename] : cluster
    -e [filename] : expression tsv file 
    --seurat-dir : seurat.R output direcory giving default cluster and expression table
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--seurat-dir', default='merged/seurat')
    parser.add_argument('-c')
    parser.add_argument('-e')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--open', action='store_true')
    parser.add_argument('-o', default='out')
    args, cmds = parser.parse_known_args()

    logger = tkutil.get_logger(os.path.basename(__file__))
    seurat_dir = pathlib.Path(args.seurat_dir)
    outdir = pathlib.Path(args.o)
    auto_open = args.open
    outdir.mkdir(exist_ok=True)
    fn_clu = args.c or (seurat_dir / 'seurat.clusters.tsv').as_posix()
    fn_count = args.e or (seurat_dir / 'seurat.count.tsv').as_posix()

    clu = pd.read_csv(fn_clu, sep='\t', index_col=0)
    data = _load_mt_counts(fn_count)
    cluster = clu.values[:,0].reshape(-1)

    n_clusters = np.max(cluster) + 1
    n_cells = data.shape[1]
    table = np.zeros((n_clusters, 3), dtype=np.float32)
    traces = []
    for i in range(n_cells):
        cn = cluster[i]
        vals = data.iloc[:,i].values.reshape(-1)
        # print(cn, vals)
        if cn >= 0:
            table[cn] += [vals[1] / (vals[0] + vals[1]), vals[0] + vals[1], 1]
    df = pd.DataFrame(table, columns=['mt_ratio', 'counts', 'n'])
    avg = df.copy()
    ncells = df.iloc[:,2].values.reshape(-1)
    cnames = []
    colors = []
    for i in range(n_clusters):
        # cnames.append('C{}_{}'.format(i + 1, ncells[i]))
        cnames.append('C{}'.format(i+1))
        avg.iloc[i,0:2] /= ncells[i]
    maxval = np.array(np.max(avg, axis=0)).reshape(-1)
    for i in range(n_clusters):
        d1 = avg.iloc[i,1] / maxval[1]
        d2 = avg.iloc[i,0] / maxval[0]
        red = int((1-d1) * 255)
        green = int((d1 + d2) * 127)
        blue = int((1-d2) * 255)
        colors.append('rgb({},{},{})'.format(red, green, blue))
    traces.append(go.Scatter(x=avg.iloc[:,0], y=avg.iloc[:,1], mode='markers+text', 
        textposition='middle center',
        text=cnames, marker=dict(size=50, symbol='circle-open', color=colors, line=dict(width=10, color=colors))))
    fn_tsv = outdir / 'mt.cnt.tsv'
    fn_pdf = outdir / 'cellstate.pdf'
    fn_html = outdir / 'cellstate.html'
    pd.DataFrame(avg, columns=['mt', 'count', 'n_cells'], index=cnames).to_csv(fn_tsv.as_posix(), sep='\t')

    fig = go.Figure(traces)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(linewidth=1, linecolor='black', mirror='ticks', title='MT ratio', range=[0,0.15])
    fig.update_yaxes(linewidth=1, linecolor='black', mirror='ticks', title='Mean read count')
    plotly.offline.plot(fig, filename=fn_html.as_posix(), auto_open=auto_open)
    fig.write_image(fn_pdf.as_posix())

def draw_umap_scatter():
    """UMAP scatter plot with cluster/study filter
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--seurat-dir', default='merged/seurat', metavar='directory', help='seurat.R script output directory')
    parser.add_argument('-c', default=None, metavar='cluster file', help='precalculated cluster fiile in tsv format')
    parser.add_argument('--cluster', nargs='+', default=['C7', 'C11', 'C1', 'C9', 'C10', 'C2'])
    parser.add_argument('-e', default=None, metavar='tsv file', help='Gene expression table (TSV format)')
    parser.add_argument('-u', default=None, metavar='tsv file', help='Coordinates file in tsv')
    parser.add_argument('--verbose', action='store_true', help='verbosicty')
    parser.add_argument('--no-pdf', action='store_true', help='PDF is not necessary')
    parser.add_argument('--open', action='store_true', help='Open HTML after drawing a chart')    
    parser.add_argument('-o', default='out', metavar='directory', help='output directory')
    parser.add_argument('--prefix', default='')
    parser.add_argument('--color-schema', default='D3', metavar='string', help='color schema by name (Vivid, Safe, Prism, ...)')
    parser.add_argument('--marker-size', type=int, default=4, metavar='integer', help='scatter plot marker size')
    args, cmds = parser.parse_known_args()
    verbose = args.verbose
    cs_ = args.color_schema
    cs_ = re.sub('\\W', '', cs_)
    color_schema = eval('px.colors.qualitative.{}'.format(cs_))

    logger = tkutil.get_logger(os.path.basename(__file__))
    if verbose: logger.setLevel(10)
    auto_open = args.open
    nopdf = args.no_pdf
    seurat_dir = pathlib.Path(args.seurat_dir)
    outdir = pathlib.Path(args.o)
    outdir.mkdir(exist_ok=True)
    fn_clu = args.c or (seurat_dir / 'seurat.clusters.tsv').as_posix()
    fn_coord = args.u or (seurat_dir / 'seurat.umap.tsv').as_posix()
    fn_chart = outdir / 'clusters.html'
    fn_pdf = outdir / 'clusters.pdf'
    prefix = args.prefix
    if prefix is not None and len(prefix) > 0:
        fn_chart = outdir / 'clusters.{}.html'.format(prefix)
        fn_pdf = outdir / 'clusters.{}.pdf'.format(prefix)

    clu = pd.read_csv(fn_clu, sep='\t', index_col=0)
    logger.info('{} clusters loaded'.format(clu.shape[0]))
    coord = pd.read_csv(fn_coord, sep='\t', index_col=0)
    logger.info('{}x{} coordinates loaded'.format(coord.shape[0], coord.shape[1]))
    marker_size = args.marker_size

    if clu.shape[0] != coord.shape[0]:
        raise Exception('barcode of cluster and coordinates are incompatible {} / {}'.format(clu.shape[0], coord.shape[0]))

    if prefix is not None:
        cell_filter = re.compile('^' + prefix)
    else:
        cell_filter = re.compile('')
    if args.cluster is None:
        clusters = list(range(0, n_clusters))
    else:
        clusters = []
        for c in args.cluster:
            m = re.match('C(\\d+)', c)
            if m:
                clusters.append(int(m.group(1))-1)
            elif c.isdigit():
                clusters.append(int(c))
            else:
                raise Exception('cluster {} is malformed cluster label'.format(c))
    cluval = clu.values[:,0].reshape(-1)
    traces = []
    traces.append(
        go.Scattergl(x=coord.iloc[:,0],y=coord.iloc[:,1], 
            name='Total {} cells'.format(coord.shape[0]), mode='markers',
            marker=dict(size=2, color='lightgray')
            )
        )

    for j, c in enumerate(clusters):
        idx = []
        for i, cn in enumerate(cluval):
            if c == cn:
                barcode = clu.index[i]
                if cell_filter.match(barcode):
                    idx.append(i)
        xy = coord.values[idx,:]
        cname = prefix + 'C{} ({})'.format(c + 1, len(idx))
        traces.append(go.Scattergl(x=xy[:,0], y=xy[:,1], name=cname, mode='markers',
            marker=dict(size=marker_size, color=color_schema[j % len(color_schema)])))
    fig = go.Figure(traces)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(linewidth=1, linecolor='black', mirror='ticks', title=coord.columns[0])
    fig.update_yaxes(linewidth=1, linecolor='black', mirror='ticks', title=coord.columns[1])
    logger.info(f'writing HTML file to {fn_chart}')
    plotly.offline.plot(fig, filename=fn_chart.as_posix(), auto_open=auto_open)
    if not nopdf:
        logger.info(f'writing PDF file to {fn_pdf}')
        fig.write_image(fn_pdf)

def __load_degs(filename):
    pass

def interpret_colorschema(colors, default_color='D3'):
    """Convert color codes to color schema.
    Plotly, Set1, D3... -> RGB arrays
    AAA,BBB,CCC -> rgb
    """
    import plotly.colors
    if colors is None:
        colors = default_color
    if colors.find(',') >= 0: # customized color
        colorschema = []
        for item in colors.split(','): # HEX color
            # print(item)
            item = item.strip('#')
            if len(item) == 3:
                item = item[0] + item[0] + item[1] + item[1] + item[2] + item[2]
            val = int(item, base=16)
            red = int(val >> 16) & 255
            green = int(val >> 8) & 255
            blue = val & 255
            colorschema.append(f'rgb({red},{green},{blue})')
    else: # color by name
        try:
            colorschema = eval(f'plotly.colors.sequential.{colors}')
        except:
            raise Exception(f'cannot interpret color name {colors}')
    return colorschema

def draw_heatmaps(argument=None, logger=None):
    """Heatmap
    --seurat-dir : directory to seurat.R output
    -c : cluster
    -e : expression
    --color : color name (plotly name)
    --zmax : z
    """
    if argument is None:
        parser = argparse.ArgumentParser()
        parser.add_argument('--seurat-dir', default='merged/seurat', help='single cell counts')
        # parser.add_argument('-i', nargs='+', default=['merged/analysis/nkt1.genes.tsv', ], help='cluster analysis file, output of get_prominent')
        # parser.add_argument('-q', default=0.01, type=float)
        parser.add_argument('-e', type=float)
        parser.add_argument('-c', type=float)
        # parser.add_argument('-r', default=2.0, type=float)
        # parser.add_argument('-p', default=0.1, type=float)
        parser.add_argument('--verbose', action='store_true')
        parser.add_argument('--unsort', action='store_true')
        # parser.add_argument('-d', default=None, help='expression by cluster file')
        parser.add_argument('-o', default='heatmaps', metavar='directory', help='output directory')
        parser.add_argument('--exclude-genes', nargs='+', default=[], help='exclude genes')
        parser.add_argument('--open', action='store_true', default=False, help='auto open heatmaps')
        parser.add_argument('-g', nargs='+', default=None)
        parser.add_argument('--color', default='RdBu_r')
        parser.add_argument('--zmax', type=float, default=4)
        parser.add_argument('--singlecell', action='store_true')
        parser.add_argument('--method', choices=['zscore', 'rawcount'], default='zscore')

        parser.add_argument('--cluster', nargs='+', default=['C7', 'C11', 'C1', 'C9', 'C10', 'C2'])

        args, cmds = parser.parse_known_args()
    else:
        args = arguments

    seurat_dir = pathlib.Path(args.seurat_dir)
    outdir = pathlib.Path(args.o)
    # phi_thr = args.p
    use_sort = not args.unsort
    zmax = args.zmax
    # qthr = args.q
    # min_expr = args.e
    # min_ratio = args.r
    scmode = args.seurat_dir is not None
    excluded_genes = args.exclude_genes
    outdir.mkdir(exist_ok=1)
    fn_clu = args.c or (seurat_dir / 'seurat.clusters.tsv').as_posix()
    fn_expr = args.e or (seurat_dir / 'seurat.normalized.tsv').as_posix()
    fn_log = outdir / 'heatmap.log'
    fn_info = outdir / 'heatmap.info'
    fn_chart = outdir / 'heatmap.html'
    fn_tsv = outdir / 'heatmap.tsv'
    fn_raw = outdir / 'heatmap.raw.tsv'
    fn_sc_tsv = outdir / 'heatmap.sc.tsv'
    method = args.method

    genes = args.g
    selected_clusters = args.cluster
    clusters = args.cluster
    # colorschema = args.color

    logger = logger or tkutil.get_logger(os.path.basename(__file__), fn_log)
    if args.verbose:
        logger.setLevel(10)

    colorschema = interpret_colorschema(args.color)
    logger.info('color ' + str(colorschema))

    # fn_data = args.d or args.i[0]
    # fn_out = args.o
    common_genes = set()
    info = {
        'command':sys.argv,
        # 'input_data':args.i,
        'expression':fn_expr,
        'cluster':fn_clu,
        'color':colorschema,
        # 'threshold':
        # {
        #     'min_expression':min_expr,
        #     'min_ratio':min_ratio,
        #     'qvalue':qthr,
        #     'phi':phi_thr,
        # }
        'zmax':zmax,
        'expression':fn_expr,
    }
    # Select genes
    # nkt1.genes.tsv
    if genes is None:
        logger.info('select genes by phi correlation')
        raise Exception('not implemented')
        for fn in args.i:
            genes = []
            with open(fn) as fi:
                items = fi.readline().strip().split('\t')
                n_clusters = len(items) - 9
                cnames = items[9:]
                for line in fi:
                    items = line.strip().split('\t')
                    if items[0].startswith('mt'): continue
                    q1 = float(items[5])
                    q2 = float(items[8])
                    e1 = float(items[2])
                    e2 = float(items[3])
                    phi = float(items[6])
                    if q1 < qthr and q2 < qthr and phi > phi_thr and e1 > min_expr and e1 > e2 * min_ratio:
                        genes.append(items[0])
            if len(common_genes) == 0:
                common_genes = set(genes)
            else:
                common_genes = set([g for g in genes if g in common_genes])
        if excluded_genes is not None:
            accepted = []
            for g in common_genes:
                if g not in excluded_genes:
                    accepted.append(g)
            common_genes = set(accepted)
        if args.genes is not None: # determine genes from analysis results
            common_genes = set(args.genes)#[g for g in args.genes if g in common_genes])
        logger.info('{} genes, {} clusters'.format(len(common_genes), n_clusters))
        genes = list(common_genes)
    else:
        # pass
        pass
    info['genes'] = genes
    logger.info('loading clusters')
    clu = pd.read_csv(fn_clu, sep='\t', index_col=0)
    cluval = clu.values.reshape(-1)
    n_clusters = np.max(cluval) + 1
    clu2idx = {}
    for i in range(clu.size):
        c = cluval[i]
        if c not in clu2idx: clu2idx[c] = []
        clu2idx[c].append(i)
    logger.info('{} clusters of {} cells detected'.format(len(clu2idx), clu.shape[0]))
    n_cells = 0
    # print(clusters)
    if clusters is None:
        for c in range(n_clusters):
            if c in clu2idx:
                n_cells += len(clu2idx[c])
        clusters = list(range(0, n_clusters))
    else:
        cl_ = []
        for c_ in clusters:
            m = re.match('C(\\d+)', c_)
            if m:
                c = int(m.group(1)) - 1
            elif c_.isdigit():
                c = int(c_)
            else:
                raise Exception(f'bad cluster name {c_}')
            if c not in clu2idx:
                raise Exception('no cluster ID={}'.format(c))
            cl_.append(c)
            n_cells += len(clu2idx[c])
        clusters = cl_
    n_clusters = len(clusters)

    logger.info('{} clusters will be used'.format(len(clusters)))

    logger.info('loading expression')
    expr = cellrangerwrapper.load_gene_expression(fn_expr, genes=genes, logger=logger)


    genes = [g for g in genes if g in expr.index]
    n_genes = len(genes)
    logger.info('{} x {} matrix loaded'.format(len(genes), expr.shape[1]))

    # sort by cluster
    avgmat = np.zeros((n_genes, n_clusters))
    scmat = np.zeros((n_genes, n_cells))

    avgcol = []
    sccol = []
    idx = []
    for ci, c in enumerate(clusters):
        # print(ci, c, clu2idx[c])
        submat = expr.iloc[:,clu2idx[c]]
        a_ = np.array(np.mean(submat, axis=1)).reshape(-1, 1)
        # print(ci, c, a_.shape)
        avgmat[:,ci:ci+1] = a_
        scmat[:,len(idx):len(idx) + submat.shape[1]] = submat
        cname = 'C{}'.format(c + 1)
        avgcol.append(cname)
        for i in range(submat.shape[1]):
            sccol.append('{}.{:04d}'.format(cname, i + 1))
        idx += clu2idx[c]
    df_sc = pd.DataFrame(scmat, index=genes, columns=sccol)

    a = np.mean(avgmat, axis=1).reshape(-1, 1)
    s = np.std(avgmat, axis=1).reshape(-1, 1)
    s[s<=0.1] = 0.1

    diffmat = avgmat - a
    zscore = diffmat / s
    zscore[np.isnan(zscore)] = 0
    df_avg = pd.DataFrame(zscore, index=genes, columns=avgcol)

    if use_sort:
        import sklearn.cluster
        clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)
        clu.fit(df_avg)
        scores = [[] for i in range(n_clusters)]
        for i in range(n_genes):
            vals = zscore[i,:]
            if np.std(vals) < 0.1: continue
            l = clu.labels_[i]
            if l >= 0:
                scores[l].append(np.argmax(vals))
        weights = [np.mean(scores[i]) for i in range(n_clusters)]
        idx = []
        for row in sorted(range(n_genes), key=lambda j:weights[clu.labels_[j]]):
            idx.append(row)
        # print(idx)
        df_avg = df_avg.iloc[idx,:]
        df_sc = df_sc.iloc[idx,:]
        df_geneclus = pd.DataFrame(clu.labels_[idx], index=df_avg.index, columns=['Cluster', ])
        df_geneclus.to_csv((outdir / 'sorted.genes.tsv').as_posix(), sep='\t')
        info['gene_cluster'] = dict(zip(df_geneclus.index, [int(x_) for x_ in clu.labels_[idx]]))
    else:
        # a = np.mean(scmat, axis=1).reshape(-1, 1)
        # s = np.std(scmat, axis=1).reshape(-1, 1)
        # s[s<=0] = 1
        # diffmat = scmat - a
        # zscore = diffmat / s
        # zscore[np.isnan(zscore)] = 0
        df_sc = pd.DataFrame(scmat, index=genes, columns=sccol)
        # print(df_sc)

    figs = plotly.subplots.make_subplots(cols=2, rows=1, subplot_titles=['cluster average', 'by cell'])
    figs.add_trace(go.Heatmap(z=df_avg.values, x=avgcol, y=df_avg.index, zmin=-zmax, zmax=zmax, colorscale=colorschema), row=1, col=1)
    figs.add_trace(go.Heatmap(z=df_sc.values, x=sccol, y=df_sc.index, zmin=-zmax, zmax=zmax, colorscale=colorschema), row=1, col=2)
    figs['layout']['yaxis1']['autorange'] = 'reversed'
    figs['layout']['yaxis2']['autorange'] = 'reversed'
    plotly.offline.plot(figs, filename=fn_chart.as_posix())
    df_avg.to_csv(fn_tsv.as_posix(), sep='\t', float_format='%.3f')
    df_sc.to_csv(fn_sc_tsv.as_posix(), sep='\t', float_format='%.3f')
    # df_raw.to_csv(fn_raw.as_posix(), sep='\t', float_format='%.3f')
    with open(fn_info, 'w') as fo:
        json.dump(info, fo, indent=2)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(__file__)
    args, cmd = parser.parse_known_args()
    if len(cmd) == 0 or cmd[0] == 'swarm':
        swarmplot_by_cluster()
    elif cmd[0] == 'heatmap':
        draw_heatmaps()
    elif cmd[0] == 'correlation':
        draw_correlation()
    elif cmd[0] == 'cellstate':
        draw_cell_state()
    elif cmd[0] == 'umap':
        draw_umap_scatter()
    else:
        raise Exception('unknown command' + cmd[0])



