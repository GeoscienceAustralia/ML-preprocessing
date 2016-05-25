"""
Run cross-validation metrics on a model prediction

.. program-output:: validatemodel --help
"""
import logging
import sys
import json
import os.path
import click as cl
import numpy as np
import matplotlib.pyplot as pl

from sklearn.metrics import r2_score, explained_variance_score
from revrand.metrics import smse, mll, msll, lins_ccc

import uncoverml.defaults as df
from uncoverml import geoio, feature
from uncoverml.validation import input_cvindex, input_targets
from uncoverml.models import apply_multiple_masked

log = logging.getLogger(__name__)


def get_first_dim(Y):

    return Y[:, 0] if Y.ndim > 1 else Y


# Decorator to deal with probabilistic output for non-probabilistic scores
def score_first_dim(func):

    def newscore(y_true, y_pred, *args, **kwargs):

        return func(y_true.flatten(), get_first_dim(y_pred), *args, **kwargs)

    return newscore

metrics = {'r2_score': r2_score,
           'expvar': explained_variance_score,
           'smse': smse,
           'lins_ccc': lins_ccc,
           'mll': mll,
           'msll': msll
           }

probscores = ['msll', 'mll']


@cl.command()
@cl.option('--quiet', is_flag=True, help="Log verbose output",
           default=df.quiet_logging)
@cl.option('--outfile', type=cl.Path(exists=False), default=None,
           help="File name (minus extension) to save output too")
@cl.option('--plotyy', is_flag=True, help="Show plot of the target vs."
           "prediction, otherwise just save")
@cl.argument('cvindex', type=(cl.Path(exists=True), int))
@cl.argument('targets', type=cl.Path(exists=True))
@cl.argument('prediction_files', type=cl.Path(exists=True), nargs=-1)
def main(cvindex, targets, prediction_files, plotyy, outfile, quiet):
    """
    Run cross-validation metrics on a model prediction.

    The following metrics are evaluated:

    - R-square
    - Explained variance
    - Standardised Mean Squared Error
    - Lin's concordance correlation coefficient
    - Mean Gaussian negative log likelihood (for probabilistic predictions)
    - Standardised mean Gaussian negative log likelihood (for probabilistic
      predictions)
    """

    # setup logging
    if quiet is True:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # Read cv index and targets
    cvind = input_cvindex(cvindex[0])
    s_ind = np.where(cvind == cvindex[1])[0]
    t_ind = np.where(cvind != cvindex[1])[0]

    Y, target_indices = input_targets(targets)
    Y = Y[target_indices]

    Yt = Y[t_ind]
    Ys = Y[s_ind]
    Ns = len(Ys)

    # build full filenames
    full_filenames = [os.path.abspath(f) for f in prediction_files]
    log.debug("Input files: {}".format(full_filenames))

    # verify the files are all present
    files_ok = geoio.file_indices_okay(full_filenames)
    if not files_ok:
        log.fatal("Input file indices invalid!")
        sys.exit(-1)

    # Load all prediction files
    filename_dict = geoio.files_by_chunk(full_filenames)
    pred_dict = feature.load_data(filename_dict, range(len(filename_dict)))

    # Deal with missing data
    EYs = feature.data_vector(pred_dict)

    # See if this data is already subset for xval
    if len(EYs) > Ns:
        EYs = EYs[s_ind]

    scores = {}
    for m in metrics:

        if m not in probscores:
            score = apply_multiple_masked(score_first_dim(metrics[m]),
                                          (Ys, EYs))
        elif EYs.ndim == 2:
            if m == 'mll' and EYs.shape[1] > 1:
                score = apply_multiple_masked(mll, (Ys, EYs[:, 0], EYs[:, 1]))
            elif m == 'msll' and EYs.shape[1] > 1:
                score = apply_multiple_masked(msll, (Ys, EYs[:, 0], EYs[:, 1]),
                                              (Yt,))
            else:
                continue
        else:
            continue

        scores[m] = score
        log.info("{} score = {}".format(m, score))

    if outfile is not None:
        with open(outfile + ".json", 'w') as f:
            json.dump(scores, f, sort_keys=True, indent=4)

    # Make figure
    if plotyy or outfile is not None:
        fig = pl.figure()
        maxy = max(Ys.max(), get_first_dim(EYs).max())
        miny = min(Ys.min(), get_first_dim(EYs).min())
        apply_multiple_masked(pl.plot, (Ys, get_first_dim(EYs)), ('k.',))
        pl.plot([miny, maxy], [miny, maxy], 'r')
        pl.grid(True)
        pl.xlabel('True targets')
        pl.ylabel('Predicted targets')
        pl.title('True vs. predicted target values.')
        if outfile is not None:
            fig.savefig(outfile + ".png")
        if plotyy:
            pl.show()