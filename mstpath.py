import os, sys, re, argparse
import tkmst
import pandas as pd
import numpy as np
import logging
import scipy.stats
import numpy as np
import gzip, hashlib
import plotly.graph_objs as go
import plotly
import plotly.io
import mygo

def __get_logger(logger=None, logfile=None):
    if logger is None:
        logger = logging.getLogger(sys._getframe().f_code.co_name)
    # set logging
    def _set_log_handler(logger, handler):#, verbose):
        handler.setFormatter(logging.Formatter('%(asctime)s %(name)s:%(lineno)s %(funcName)s [%(levelname)s]: %(message)s'))
        logger.addHandler(handler)
        return logger
    _set_log_handler(logger, logging.StreamHandler())
    if logfile is not None:
        _set_log_handler(logger, logging.FileHandler(logfile))
    # logger.setLevel(logging.ERROR)
    logger.propagate = False
    return logger

class MSTNode(object):
    def __init__(self, name, size=100, layer=0, number=0):
        self.name = name
        self.size = size
        self.layer = layer
        self.number = number

def draw_mst_path(nodes, paths, outdir, **kwargs):
    import reportlab.pdfgen.canvas
    # xstep = kwargs.get('xstep', 50)
    ystep = kwargs.get('ystep', 30)
    xpos = kwargs.get('x', 400)
    ypos = kwargs.get('y', 400)
    circle_size = kwargs.get('size', 0.4)
    labels = kwargs.get('labels', None)
    fn_fig = os.path.join(outdir, 'chart.pdf')
    fn_tree = os.path.join(outdir, 'path.tsv')
    cnv = reportlab.pdfgen.canvas.Canvas(fn_fig)

    coordinates = []
    # curves = []    
    layers = {}
    for i, node in enumerate(nodes):
        if node.layer not in layers:
            layers[node.layer] = 0
        layers[node.layer] += 1
    max_size = max(layers.values())

    layer_count = {}
    for l in layers.keys(): layer_count[l] = 0
    xstep = kwargs.get('xstep', 550 / len(layers))
    for i, node in enumerate(nodes):
        layer_size = layers[node.layer]
        y = ypos - (layer_count[node.layer] - layer_size / 2) * ystep 
        layer_count[node.layer] += 1
        x = (node.layer - len(layers) / 2 - 1) * xstep + xpos
        r = np.sqrt(node.size + 25) * circle_size
        coordinates.append((x, y, r))

    cnv.setFont('Helvetica', 10)
    if labels is not None:
        for i, l in enumerate(labels):
            x = (i - len(layers) / 2 - 1) * xstep + xpos
            y = ypos - max_size * 0.5 * ystep + 10
            cnv.drawString(x - cnv.stringWidth(l) / 2, y, l)
    cnv.setStrokeColor('gray')
    cnv.setFont('Helvetica', 6)
    paramdist = []
    for path in paths:
        if isinstance(path, list) or isinstance(path, tuple) and len(path) > 2:
            d = path[2]
            paramdist.append(d)
    #     else:
    #         print(path)
    # print(paths)
    # print(paramdist)
    # exit()
    dmin = min(paramdist)
    dmax = max(paramdist)

    for i, path in enumerate(paths):
    # for i, node in enumerate(nodes):
        pathinfo = paths[i]
        if isinstance(pathinfo, int):
            node = i
            parent = pathinfo
            param = 0
        else:
            node = pathinfo[0]
            parent = pathinfo[1]
            if len(pathinfo) > 2:
                param = pathinfo[2]
        # parent = paths[i]
        if parent >= 0:
            cnv.setLineWidth(0.5)
            x0, y0 = coordinates[node][0:2]
            x1, y1 = coordinates[parent][0:2]
            if param >= 0:
                degree = (param - dmin) / (dmax - dmin)
                # cnv.setLineWidth(0.5)#max(0.5, 2 - param))
                red = 1.0
                green = blue = 0.5
                red = min(1.0, max(0, 0.5 + (0.5 - param * 10)))
                # green = blue = max(0, param)
                cnv.setLineWidth(2.5 - degree * 2)
                cnv.setStrokeColorRGB(red, green, blue)
                print(i, parent, param, red, green, blue)

            # curves((x0, y0, x1, y1))
            if x0 == x1: 
                dist = y1 - y0
                if (dist > 0) == (i % 2 == 0):
                    x2 = x3 = x0 - dist / 10
                else:
                    x2 = x3 = x0 + dist / 10
                y2 = y0 + dist / 4
                y3 = y1 - dist / 4
                # cnv.setLineWidth(0.25)
                cnv.bezier(x0, y0, x2, y2, x3, y3, x1, y1)
            else:
                # cnv.setLineWidth(0.6)
                cnv.line(x0, y0, x1, y1)
    cnv.setLineWidth(1)
    for i, node in enumerate(nodes):
        x, y, r = coordinates[i]
        cnv.setStrokeColor('black')
        cnv.setFillColor('white')
        cnv.circle(x, y, r, 1, 1)
        cnv.setFillColor('black')
        # print(x, y, r, node.name)
        cnv.drawString(x - cnv.stringWidth(node.name) / 2, y - 3, node.name)
        sstr = str(node.size)
        cnv.drawString(x - cnv.stringWidth(sstr) / 2, y - 10, sstr)
    cnv.save()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seurat-dir', default=None)
    parser.add_argument('-c' ,default=None)
    parser.add_argument('-i', default=None, help='PCA or UMAP file, row index is barcode')
    parser.add_argument('-u', default=None, help='Graph coordinates')
    # parser.add_argument('-i', nargs='+', help='anova.tsv files calculated by seurat.R')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--forced', action='store_true')
    parser.add_argument('--metrics', choices=['pearson', 'cosine'], default='pearson')
    parser.add_argument('-n', default=50, type=int, help='degree of freedom')
    parser.add_argument('--min-cluster-size', type=int, default=0)
    parser.add_argument('--cluster-filter', type=int, default=None, nargs='+', 
        help='designate analyzing clusters by index')
    parser.add_argument('--organism', default='human')
    parser.add_argument('--geneset', default=None)
    parser.add_argument('--open', action='store_true')
    # parser.add_argument('--alias')

    parser.add_argument('-o', default='mst')
    args = parser.parse_args()
    auto_open = args.open
    dof = args.n
    minimum_cluster_size = args.min_cluster_size
    cluster_filter = args.cluster_filter

    outdir = args.o
    forced = args.forced
    os.makedirs(outdir, exist_ok=True)
    genes = None

    logger = __get_logger(logging.getLogger(sys._getframe().f_code.co_name))
    if args.verbose: logger.setLevel(logging.DEBUG)

    metrics = args.metrics
    if metrics == 'pearson':
        __metrics = lambda x, y: 1 - scipy.stats.pearsonr(x, y)[0]
    elif metrics == 'cosine':
        __metrics = lambda x, y: 1 - np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y)
    else:
        raise Exception('unknown metrics {}'.format(metics))

    basedir = args.seurat_dir
    if basedir is not None:
        fn_cluster = os.path.join(basedir, 'seurat.clusters.tsv')
        fn_data = os.path.join(basedir, 'seurat.pca.tsv')
        fn_coord = os.path.join(basedir, 'seurat.umap.tsv')
    if args.c:
        fn_cluster = args.c
    if args.i:
        fn_data = args.i
    if args.u:
        fn_coord = args.u

    if args.geneset is not None:
        if os.path.exists(args.geneset):
            with open(args.geneset) as fi:
                for line in fi:
                    if line.startswith('#'): continue
                    genes.append(line.strip().split('\t')[0])
            genes = sorted(genes)
        else:
            organism = args.organism
            genes = None
            for item in args.geneset.split('+'):
                logger.info('GO {} : {}'.format(organism, item))

                genes_ = mygo.get_genes_having_go(organism=organism, keyword=item)
                if genes is None:
                    genes = set(genes_)
                else:
                    genes = set([g for g in genes_ if g in genes])
            genes = list(sorted(genes))
            logger.info('{} genes were selected'.format(len(genes)))


    md5 = hashlib.md5()
    # alias = {}
    md5.update(metrics.encode('utf-8'))
    # for f in args.i:
    md5.update('cluster={},data={},coord={},dof={},metrics={};min_size={};'.format(
        os.path.abspath(fn_cluster), 
        os.path.abspath(fn_data), 
        os.path.abspath(fn_coord), dof, metrics,
        minimum_cluster_size)
        .encode('utf-8'))
    if cluster_filter is not None:
        md5.update('clusterfilter={};'.format(str(cluster_filter)).encode('utf-8'))

    if genes is not None:
        md5.update('genes={}'.format(',', list(genes)).encode('utf-8'))

    fn_dist = os.path.join(outdir, '.' + md5.hexdigest()[0:8] + '.mstmat')
    df_cls = pd.read_csv(fn_cluster, sep='\t', index_col=0)
    df_cls.columns = ['cluster', ]
    df_coord = pd.read_csv(fn_coord, sep='\t', index_col=0)
    logger.info('cluster file, {}, coord file {}'.format(str(df_cls.shape), str(df_coord.shape)))

    # filter clusters by cell number
    if cluster_filter is None and minimum_cluster_size > 0:
        cc = {}
        for i, c in enumerate(df_cls['cluster'].values.reshape(-1)):
            if c not in cc: cc[c] = []
            cc[c].append(i)
        # copied = df_cls.copy()
        ignoring_index = []
        for c in cc.keys():
            if len(cc[c]) < minimum_cluster_size:
                ignoring_index += cc[c]
                pass
        logger.info('{} cells were ignored by cluster size'.format(len(ignoring_index)))
        df_cls.iloc[ignoring_index] = -1
    
    # filter clusters by arguments
    if cluster_filter is not None:
        allowed = set(cluster_filter)
        ignoring = []
        for i, c in enumerate(df_cls['cluster'].values.reshape(-1)):
            if c not in allowed:
                ignoring.append(i)
        logger.info('{} cells were ignored because of cluter-filter'.format(len(ignoring)))
        df_cls.iloc[ignoring] = -1
        

    if os.path.exists(fn_dist) and not forced:
        distance_matrix = pd.read_csv(fn_dist, sep='\t', index_col=0)
        n_clusters = distance_matrix.shape[0]
    else:
        df_data = pd.read_csv(fn_data, sep='\t', index_col=0)
        # if genes is not None:
        #     available_genes = [g for g in df_data.index if g in genes]
        #     df_data = df_data.loc[available_genes]
        #     logger.info('using {}x{} matrix'.format(df_data.shape[0], df_data.shape[1]))

        n_clusters = np.max(df_cls.values) + 1
        logger.info('{} clusters'.format(n_clusters))
        c2pos = []
        if df_data.shape[1] > dof:
            df_data = df_data.iloc[:,0:dof]

        dmat = np.zeros((n_clusters, n_clusters), dtype=np.float32)
        cluster_names = []
        for c in range(n_clusters):
            barcodes = df_cls[df_cls['cluster']==c].index
            if len(barcodes) > 0:
                pos = np.mean(df_data.loc[barcodes], axis=0)
            else:
                pos = None
            cluster_names.append('C{}_{}'.format(c + 1, len(barcodes)))
            if pos is not None:
                logger.info('Distance of {}'.format(cluster_names[-1]))
                # if len(barcodes) > 0:
                j = len(c2pos)
                for i, p in enumerate(c2pos):
                    if p is None: continue
                    if len(pos) > 0 and len(p) > 0:
                        dist = __metrics(pos, p)
                        dmat[i,j] = dmat[j,i] = dist
                    else:
                        dmat[i,j] = dmat[j,i] = 0
            c2pos.append(pos)
            # c2pos.append(pos)
        distance_matrix = pd.DataFrame(dmat, columns=cluster_names, index=cluster_names)
        distance_matrix.to_csv(fn_dist, sep='\t', float_format='%.6f')
    
    # skip clusters without cells
    available = []
    for c in range(n_clusters):
        if df_cls[df_cls['cluster']==c].shape[0] == 0: # no cells
            logger.info('skip cluster {}'.format(c))
        else:
            available.append(c)
    # print(distance_matrix)
    distance_matrix = distance_matrix.iloc[available, available]
    # print(distance_matrix)
    # exit()

    # points = []
    positions = {}
    traces = []    
    # unclustered
    barcodes = df_cls[df_cls['cluster']==-1].index
    if len(barcodes) > 0:
        df_xy = df_coord.loc[barcodes]
        print(df_xy.shape)
        print(df_xy.values[:,0])
        print(df_xy.values[:,1])
        traces.append(
            go.Scattergl(x=df_xy.values[:,0], y=df_xy.values[:,1],
                mode='markers', name='others' ,marker=dict(size=2, color='lightgray'))
            )
    # clusters
    cn2idx = {}
    for c in range(n_clusters):
        barcodes = df_cls[df_cls['cluster']==c].index
        if len(barcodes) == 0: continue
        cn2idx[len(cn2idx)] = c
        df_xy = df_coord.loc[barcodes]
        pos = np.mean(df_xy, axis=0)
        positions[c] = pos#.append(pos)
        x = df_xy.iloc[:,0].values
        y = df_xy.iloc[:,1].values
        traces.append(go.Scattergl(x=x, y=y, name='C{}'.format(c + 1), mode='markers', marker=dict(size=3)))

    paths = []
    dmat_ = distance_matrix.values
    print(distance_matrix)
    print(positions)
    for i, upper in enumerate(tkmst.find_mst_m(distance_matrix.values)):
        c1 = cn2idx.get(i, -1)
        c2 = cn2idx.get(upper, -1)
        if upper < 0:
            paths.append((i, -1, 0))
        else:
            logger.info('{}=>{}:{}\t{}=>{}:{}\t{}'.format(
                i, c1, distance_matrix.columns[i], 
                upper, c2, distance_matrix.columns[upper], dmat_[upper, i]))
            paths.append((i, upper, dmat_[upper, i]))
            x0, y0 = positions[c1][0:2]
            x1, y1 = positions[c2][0:2]
            logger.info('({:.0f},{:.0f}) - ({:.0f},{:.0f})'.format(x0, y0, x1, y1))

            traces.append(
                go.Scattergl(x=[x1, x0], y=[y1, y0], mode='lines+markers+text', 
                text=['C{}'.format(c2+ 1), 'C{}'.format(c1+1)], line=dict(color='black', width=3),
                marker=dict(size=32, symbol='circle', color='white', line=dict(color='black', width=2)))
            )

    fig = go.Figure(traces)
    fig.update_layout(plot_bgcolor='white', showlegend=False)
    fig.update_xaxes(linecolor='black', linewidth=1, mirror='ticks', title='UMAP_1')
    fig.update_yaxes(linecolor='black', linewidth=1, mirror='ticks', title='UMAP_2')

    fn_html = os.path.join(outdir, 'scattermap.html')
    fn_pdf = os.path.join(outdir, 'scattermap.pdf')
    plotly.offline.plot(fig, filename=fn_html, auto_open=auto_open)
    fig.update_layout(width=800, height=800)
    plotly.io.write_image(fig, fn_pdf, format='pdf')

    exit()


    nodes = []
    import collections
    layers = collections.OrderedDict()
    # layers = {}
    for c in distance_matrix.columns:
        m = re.match('(\\w+)\\.(C\\d+)_(\\d+)', c)
        if m:
            size = int(m.group(3))
            cluster = m.group(2)
            sample = m.group(1)
            l = sample
        else:
            l = c.split('.')[0]
            size = 100
            cluster = l
        if l not in layers: layers[l] = len(layers)
        n = MSTNode(cluster, size, layers[l])
        # print(c, cluster, size, layers[l])
        nodes.append(n)
    # exit()
    paths = []
    dmat_ = distance_matrix.values
    for i, upper in enumerate(tkmst.find_mst_m(distance_matrix.values)):
        if upper < 0:
            paths.append((i, -1, 0))
        else:
            paths.append((i, upper, dmat_[upper, i]))
    connections = []    

    for i in range(dmat_.shape[0]):
        for j in range(i):
            # if dmat_[i,j] != dmat_[j,i]:
            d = max(dmat_[i,j], dmat_[j, i])
            print(i, j, d)
            if d < 1:#0.05:
                connections.append((i, j, d))
    #             print(i, j, dmat_[i,j], dmat_[j,i])
    # # exit()
    # # print(paths)
    # # print(connections)
    # # exit()
    # print(connections)
    # exit()

    # print(paths)
    cluster_names = distance_matrix.columns
    with open(os.path.join(outdir, 'path.tsv'), 'w') as fo:
        for i, p in enumerate(paths):
            fo.write('{}\t{}\t{}\t{}\n'.format(i, cluster_names[i], p[1], cluster_names[p[1]]))
    draw_mst_path(nodes, connections, outdir, labels=[re.sub('^\\d+', '', l) for l in layers.keys()])
    pass


if __name__ == '__main__':
    main()