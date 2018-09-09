from __future__ import division

import glob
import logging
import os
import shutil

from amuse.community import units

import spdummy
import ncmod
import sputils

# model type names
dummy_type = "dummy"
oifs_type = "oifs"
dales_type = "dales"
ncbased_type = "ncfile"
dummy_gcm_type = "dummy_gcm"
dummy_les_type = "dummy_les"
ncfile_gcm_type = "spifsnc_gcm"
ncfile_les_type = "spifsnc_les"

# Logger
log = logging.getLogger(__name__)


# Factory method for model creation
# TODO: use kwargs...   
def create_model(model_type, inputdir, workdir, nprocs=1, redirect="file", channel_type="mpi", restart=False,
                 starttime=None, index=-1, qt_forcing="sp", trestart=1000000 | units.s):
    ofile = os.path.join(workdir, model_type + ".out")
    efile = os.path.join(workdir, model_type + ".err")
    if model_type == oifs_type:
        from omuse.community.oifs.interface import OpenIFS

        if not restart:
            # some files cannot be linked:
            #   fort.4 is patched before the run
            #   ncf927 is written to by openIFS. We copy it, if it exists in the input dir
            files = [f for f in glob.glob(os.path.join(inputdir, "*"))
                     if not (f.endswith("fort.4") or f.endswith("ncf927")) ]
            log.info("Linking OpenIFS input files from %s...\n Files: %s" % (inputdir, ' '.join(files)))
            log.info("Workdir:%s" % workdir)
            sputils.link_dir(files, workdir)
            shutil.copy(os.path.join(inputdir, "fort.4"), workdir)
            try:
                shutil.copy(os.path.join(inputdir, "ncf927"), workdir)
            except:
                pass

        oifs = OpenIFS(number_of_workers=nprocs,
                       redirection=redirect,
                       redirect_stdout_file=ofile,
                       redirect_stderr_file=efile,
                       channel_type=channel_type,
                       workdir=workdir)  # , debugger='gdb')

        setattr(oifs, "workdir", workdir)

        return oifs

    elif model_type == dales_type:
        from omuse.community.dales.interface import Dales

        qt_forcings = {"sp": Dales.QT_FORCING_GLOBAL,
                       "variance": Dales.QT_FORCING_VARIANCE,
                       "local": Dales.QT_FORCING_LOCAL}

        if not restart:
            files = glob.glob(os.path.join(inputdir, '*'))
            log.info("Linking Dales input files from %s..." % inputdir)
            sputils.link_dir(files, workdir)

        dales = Dales(number_of_workers=nprocs,
                      redirection=redirect,
                      redirect_stdout_file=ofile,
                      redirect_stderr_file=efile,
                      channel_type=channel_type,
                      workdir=workdir)  # , debugger='gdb')
        dales.parameters.restart_flag = restart
        dales.parameters.trestart = trestart
        if starttime is not None:
            dales.parameters.starttime = starttime
        dales.parameters.qt_forcing = qt_forcings[qt_forcing]
        setattr(dales, "workdir", workdir)
        return dales

    elif model_type == dummy_gcm_type:
        dummy = spdummy.dummy_gcm(nprocs)
        setattr(dummy, "workdir", workdir)
        return dummy

    elif model_type == dummy_les_type:
        dummy = spdummy.dummy_les(nprocs)
        setattr(dummy, "workdir", workdir)
        return dummy

    elif model_type == ncfile_gcm_type:
        return ncmod.netcdf_gcm(os.path.join(inputdir, "spifs.nc"))

    elif model_type == ncfile_les_type:
        return ncmod.netcdf_les(os.path.join(inputdir, "spifs.nc"), index)

    else:
        raise Exception("Unsupported model type: " + model_type)
