# Superparameteriation coupling code for OpenIFS <--> Dales
#
# Fredrik Jansson, Gijs van den Oord 
# 2017-2019
# 


from __future__ import division
from __future__ import print_function

import json
import logging
import os
import shutil
import threading
from Queue import Queue  # note named queue in python 3

import datetime
import numpy
import sys
import time

import modfac
import spcpl
import sputils
import spio
import spmpi
import psutil

from amuse.community import *
from amuse.rfi import channel  # to query MPI threading support

#from amuse.rfi.channel import AsyncRequestsPool
from amuse.rfi.async_request import AsyncRequestsPool

# Logger
log = logging.getLogger(__name__)

# Module configuration variables
gcm_type = "oifs"
gcm_steps = 10  # number of gcm time steps to perform
gcm_exp_name = "TEST"  # openifs experiment name
gcm_input_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../oifs-input")  # openifs input directory
gcm_run_dir = "oifs-work"  # openifs run directory
gcm_num_procs = 1  # openifs MPI tasks
gcm_redirect = "file"  # redirection for gcm
gcm_forcing_factor = 1  # scale factor for forcings upon openifs
les_type = "dales"
les_dt = 60  # les time step (<0: adaptive)
les_spinup = 0
les_spinup_steps = 1
les_spinup_forcing_factor = 1.
les_exp_name = "test"  # les experiment name
les_input_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../dales-input")  # les input directory
les_run_dir = "dales-work"  # les run directory
les_num_procs = 1  # MPI tasks per les instance
les_redirect = "file"  # redirection for les
les_forcing_factor = 1  # scale factor for forcings upon les
les_queue_threads = sys.maxint  # les run scheduling (1: all serial, > 1: nr. of concurrent worker threads)
max_num_les = -1  # Maximal number of LES instances
init_les_state = True  # initialize les instances to the openifs column state
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../spifs-output")  # Output folder
output_name = "spifs.nc"  # output netcdf file name
channel_type = "sockets"  # amuse communication type (choose from ["sockets","mpi"])
dryrun = False  # if true, only start the GCM to examine the grid.
async_evolve = True  # time step LES instances using asynchronous amuse calls instead of Python threads (experimental)
restart = False  # restart an old run
cplsurf = False  # couple surface fields
firststep = True # flag for this being the first step - experimentally used in restarts to not log the weird first step

qt_forcing = "sp"

# Model instances:
# Global circulation model
gcm_model = None
# Local large eddy simulation models
les_models = []

output_column_indices = []
output_columns = []  # tuple (index, lat, lon)

errorFlag = False  # flag raised when a worker thread generates an exception


# Writes gridpoint file for the input geometry
def save_dryrun_info(lons, lats):
    log.info("Dry run - saving grid point coordinates in gridpoints.txt.")
    points = numpy.column_stack((lons, lats))
    numpy.savetxt('gridpoints.txt', points, fmt='%10.6f')
    log.info("Dry run finished - will exit now.")
    finalize()
    sys.exit()


# Initializes the system
def initialize(config, geometries, output_geometries=None):
    global gcm_model, les_models, output_name, async_evolve, output_column_indices, output_columns
    read_config(config)

    if not restart and os.path.exists(output_dir):
      raise Exception("output dir %s exists"%output_dir) 


    # if output_name is a relative path (as by default),
    # treat it as relative to output_dir
    if not os.path.isabs(output_name):
        output_name = os.path.join(output_dir, output_name)

    # TODO: validate input parameters
    run_dir = os.path.join(output_dir, gcm_run_dir)
    if channel_type == "nospawn":
        spmpi.send_model_colors(gcm_num_procs, les_num_procs, max_num_les)
        # TODO: Replace Dales and openifs channel factory methods...
    gcm_model = gcm_init(gcm_type, gcm_input_dir, run_dir, couple_surface=cplsurf)
    les_models = []
    lons = gcm_model.longitudes.value_in(units.deg)
    lats = gcm_model.latitudes.value_in(units.deg)
    grid_indices = sputils.get_mask_indices(zip(lons, lats), geometries, max_num_les)
    output_geoms = [] if output_geometries is None else output_geometries
    output_column_indices = sputils.get_mask_indices(zip(lons, lats), output_geoms)
    
    # exclude columns with embedded LES from the output_column_indices 
    output_column_indices = list(set(output_column_indices) - set(grid_indices))
    output_columns = [(i, lats[i], lons[i]) for i in output_column_indices]

    log.info("Creating LES models in grid columns ")  # % str(grid_indices))

    for i in grid_indices:
        log.info("%7d x=%8.3f y=%8.3f" % (i, lons[i], lats[i]))

    log.info("Extra netCDF output for grid columns ")  # % str(grid_indices))
    for c in output_columns:
        log.info("%7d x=%8.3f y=%8.3f" % (c[0], c[2], c[1]))

    if dryrun:
        save_dryrun_info(lons, lats)

    local_les_input_dir = os.path.join(output_dir, 'les-input')
    if not restart:
        # Copy les input directory into run directory. Pass the local copy to the les model init.
        shutil.copytree(les_input_dir, local_les_input_dir)

    startdate = gcm_model.get_start_datetime() - datetime.timedelta(seconds=les_spinup)

    for i in grid_indices:
        instance_run_dir = os.path.join(output_dir, les_run_dir + '-' + str(i))
        les = les_init(les_type, local_les_input_dir, instance_run_dir, startdate, i)
        gcm_model.set_mask(i)  # tell GCM that a LES instance is present at this point
        les.grid_index = i
        les.lat, les.lon = lats[i], lons[i]
        les.zh_cache = les.get_zh()
        les.zf_cache = les.get_zf()
        les_models.append(les)
    
    spio.init_netcdf(output_name, gcm_model, les_models, startdate, output_columns, append=restart,
                     with_surf_vars=cplsurf)
    log.info("Successfully initialized GCM and %d LES instances" % len(les_models))

    # Switch off async in case any model doesn't support it
    async_evolve = async_evolve and reduce(lambda p, q: p and q,
                                           [getattr(m, "support_async", True) for m in [gcm_model] + les_models])

    if channel_type != "sockets":

        # the actual thread level provided by the MPI library
        log.info(
            "MpiChannel.is_multithreading_supported(): %s" % (str(channel.MpiChannel.is_multithreading_supported())))

        if not channel.MpiChannel.is_multithreading_supported():
            if not async_evolve and les_queue_threads > 1:
                log.info(
                    "Options are set to run Dales instances from separate python threads but the MPI in use does not "
                    "support multithreading. Exit.")
                sys.exit()


    gcm_model.first_half_step_done = False
    if not restart:
        numpy.random.seed(42)  # seed generator the same way every time - for repeatable simulation

        # do first half of first time step in openIFS now, so that U,V,T get initialized

        log.info("gcm.evolve_model_until_cloud_scheme() - first step")
        gcm_model.evolve_model_until_cloud_scheme()
        log.info("gcm.evolve_model_cloud_scheme() - first step")
        gcm_model.evolve_model_cloud_scheme()
        gcm_model.first_half_step_done = True  # set flag here, to avoid repeating the half step

        # spio.update_time(gcm_model.get_model_time())
        spinup_delta_t =  les_spinup / les_spinup_steps
        spio.update_time(spinup_delta_t | units.s) # time stamps are for the LES time at the end of the step

        if init_les_state:
            spcpl.gather_gcm_data(gcm_model, les_models, True)
            # Note: no surface fluxes can be fetched after a full time step
            # But now we experimentally stay in the middle of the first step,
            # and then we can fetch them

            # get the state and apply it on les as initial state
            for les in les_models:
                u, v, thl, qt, ps, ql = spcpl.convert_profiles(les)
                spcpl.set_les_state(les, u, v, thl, qt, ps)
        
            if les_spinup > 0:
              run_spinup(les_models, gcm_model, les_spinup, les_spinup_steps)

    else:  # we're doing a restart
        pass
 
    return gcm_model, les_models

# Run loop: executes nsteps time steps of the super-parametrized GCM
def run(nsteps):
    current_process = psutil.Process(os.getpid())  # get current process, for resource usage measurement
    have_work_queue = 1 < les_queue_threads < len(les_models)
    # TODO: Check whether the gcm supports another nsteps steps
    work_queue, worker_threads = None, []
    if have_work_queue:
        work_queue, worker_threads = start_worker_threads(les_queue_threads)
    # timestep models together
    for s in range(nsteps):
        step(work_queue)
        log.info('python master usage: %s' % str(current_process.memory_full_info()))
        log.info('System total: %s' % str(psutil.virtual_memory()))
        log.info('  ---- Time step done ---')
    if have_work_queue:
        stop_worker_threads(work_queue, worker_threads)


# Spinup loop: executes nsteps time steps of the super-parametrized GCM
def run_spinup(les_list, gcm, spinup_length, spinup_steps=1):
    have_work_queue = 1 < les_queue_threads < len(les_list)
    # TODO: Check whether the gcm supports another nsteps steps
    work_queue, worker_threads = None, []
    if have_work_queue:
        work_queue, worker_threads = start_worker_threads(les_queue_threads)

    iteration_length = spinup_length / spinup_steps
    for s in range(spinup_steps):
        if s == spinup_steps - 1:
            iteration_length = spinup_length - (spinup_steps - 1) * iteration_length
        step_spinup(les_list, work_queue, gcm, spinup_length=iteration_length)

    log.info('System total: %s' % str(psutil.virtual_memory()))
    log.info('  ---- Spinup done ---')
    if have_work_queue:
        stop_worker_threads(work_queue, worker_threads)


timing_file = None

def open_timing_file():
    global timing_file
    timing_file = open(output_dir + '/timing.txt', 'a')
    # at a fresh start, write a list of the LES grid points to the file
    if not restart:
        s = '# LES grid points\n'
        s += ' '.join([str(les.grid_index) for les in les_models])
        s += '\n# timing data\n'
        timing_file.write(s)


# do one gcm time step
# step les until it catches up
def step(work_queue=None):
    global timing_file,firststep,profiles
    if not timing_file:
       open_timing_file()

    # don't write to spifs.nc at the first step of a restarted run.
    # the first step seems to repeat the last step of the previous run.
    writeCDF = (not(restart and firststep))

    if firststep:
        profile={}
    t = gcm_model.get_model_time()
    delta_t = gcm_model.get_timestep()
    log.info("gcm time at start of timestep is %s" % str(t))
    # want this message before the time stepping
    # until_cloud_scheme and cloud_scheme below do not change the model time
    
    starttime = time.time()
    gcm_walltime1 = -time.time()

    if writeCDF and not firststep: 
        spio.update_time(gcm_model.get_model_time() + (les_spinup | units.s) + delta_t)
    
    try:
        if gcm_model.first_half_step_done:
            # if we already did the first half step as part of the initialization,
            # don't repeat it now.
            gcm_model.first_half_step_done = False
        else:
            log.info("gcm.evolve_model_until_cloud_scheme()")  
            gcm_model.evolve_model_until_cloud_scheme()
            log.info("gcm.evolve_model_cloud_scheme()")
            gcm_model.evolve_model_cloud_scheme()  # note: overwrites set tendencies
    except Exception as e:
        log.error("Exception when time-stepping openIFS: %s Exiting." % e.message)
        log.error(sys.exc_info())
        finalize()
        sys.exit(1)

    gcm_walltime1 += time.time()
    gcm_model.step += 1
    delta_t=gcm_model.get_timestep() #probably don't have to calculate it one more time but this is the next gcm step...

    
    gather_gcm_data_walltime = -time.time()
    spcpl.gather_gcm_data(gcm_model, les_models, cplsurf, output_column_indices, write=writeCDF)
    gather_gcm_data_walltime += time.time()

    set_les_forcings_walltime = -time.time()
    #pool = AsyncRequestsPool()
    for les in les_models:
        if not firststep:
            profile=profiles[les]
        req=spcpl.set_les_forcings(les, gcm_model,True, firststep, profile, dt_gcm=delta_t, factor=les_forcing_factor, couple_surface=cplsurf, qt_forcing=qt_forcing, write=writeCDF)
        #for r in req.values():
        #    pool.add_request(r)
    #pool.waitall()
    set_les_forcings_walltime += time.time()
    # step les models to the end time of the current GCM step = t + delta_t
    les_wall_times, profiles = step_les_models(t + delta_t, work_queue, offset=les_spinup)
    # get les state - for forcing on OpenIFS and les stats
    set_gcm_tendencies_walltime = -time.time()
    for les in les_models:
        profile=profiles[les]
        spcpl.set_gcm_tendencies(gcm_model, les, profile=profiles[les],dt_gcm=delta_t, factor=gcm_forcing_factor, write=writeCDF)
    set_gcm_tendencies_walltime += time.time()
    gcm_walltime2 = -time.time()
    gcm_model.evolve_model_from_cloud_scheme()
    gcm_walltime2 += time.time()

    
    log.info("gcm evolved to %s" % str(gcm_model.get_model_time()))
    s = ('%10.2f %6.2f %6.2f %6.2f %6.2f %6.2f' % (starttime, gcm_walltime1, gather_gcm_data_walltime, set_les_forcings_walltime, set_gcm_tendencies_walltime, gcm_walltime2)
         + ' ' + ' '.join(['%6.2f' % t for t in les_wall_times]) + '\n')
    timing_file.write(s)
    timing_file.flush()



    # sync spifs.nc now, if we have no LES models.
    # if we do have LES models, we sync elsewhere, while the LES models are busy.
    if len(les_models) == 0:
        spio.sync_root()

    firststep = False

# Initialization function
def step_spinup(les_list, work_queue, gcm, spinup_length):
    global timing_file, firststep,profiles

    if not any(les_list): return

    if not timing_file:
        open_timing_file()
    
    if firststep:
        profile={}

    if not firststep:
        # in the very first step, this has already been done in the initialization
        spio.update_time(les_list[0].get_model_time() + (spinup_length | units.s))
    
    starttime = time.time()
    
    t_les = les_list[0].get_model_time()

    set_les_forcings_walltime = -time.time()
    pool = AsyncRequestsPool()
    for les in les_list:
	if not firststep:
	    profile=profiles[les]
        req=spcpl.set_les_forcings(les, gcm, True,firststep, profile, dt_gcm=spinup_length| units.s, factor=les_spinup_forcing_factor, couple_surface=cplsurf, qt_forcing=qt_forcing)
        for r in req.values():
            pool.add_request(r)
    pool.waitall()
    set_les_forcings_walltime += time.time()
        
    # step les models
    les_wall_times , profiles = step_les_models(t_les + (spinup_length | units.s), work_queue, offset=0)
    set_gcm_tendencies_walltime = -time.time() # assign the profile writing time to the same slot as setting gcm tendencies
    for les in les_list:
        profile=profiles[les]
        spcpl.write_les_profiles(les)
    set_gcm_tendencies_walltime += time.time()

    firststep = False   

    
    gcm_walltime1 = 0
    gcm_walltime2 = 0
    gather_gcm_data_walltime = 0
    s = ('%10.2f %6.2f %6.2f %6.2f %6.2f %6.2f' % (starttime, gcm_walltime1, gather_gcm_data_walltime, set_les_forcings_walltime, set_gcm_tendencies_walltime, gcm_walltime2)
         + ' ' + ' '.join(['%6.2f' % t for t in les_wall_times]) + '\n')
    timing_file.write(s)
    timing_file.flush()


# Function for stopping gcm and all les instances
# this is called both at a normal exit and when an exception
# is generated in one of the worker threads.
# The goal then is to 1) make all threads quit, so that the job ends instead of hanging
#                     2) make all worker_threads quit as nicely as possible so that results are saved
def finalize(save_restart=True):
    if save_restart:
        log.info("Asking Dales to save restart files.")
        # save LES restart files
        for les in les_models:
            les.write_restart()

    log.info("spifs cleanup...")
    log.info("Stopping gcm...")
    try:
        gcm_model.cleanup_code()
        gcm_model.stop()
    except Exception as e:
        log.error("Exception while stopping gcm: %s" % e.message)
    log.info("Stopping LES instances...")
    for les in les_models:
        try:
            les.cleanup_code()
            les.stop()
        except Exception as e:
            log.error("Exception while stopping LES at index %d: %s" % (les.grid_index, e.message))
    spio.cdf_root.close()
    log.info("spifs cleanup done")


# Reads input parameters from input file or dictionary
def read_config(config):
    userconf = {}
    if isinstance(config, str):
        if os.path.isfile(config):
            with open(config) as f:
                userconf = json.load(f)
        else:
            log.error("Could not find input configuration file %s" % config)
    elif isinstance(config, dict):
        userconf = config
    elif userconf:
        log.error("Could not read configurations from object of type %s" % type(config))
    if userconf:
        log.info("Configuration options read: %s" % str(userconf))
    for key in userconf:
        if key in globals():
            if hasattr(userconf[key], "__call__"):
                log.info("Skipping setting module function %s to value %s..." % (key, str(userconf[key])))
                continue
            log.info("Setting module variable %s to value %s..." % (key, str(userconf[key])))
            globals()[key] = userconf[key]


# Creates and initialized the GCM
def gcm_init(gcmtype, inputdir, workdir, couple_surface):
    typekey = gcmtype
    if gcmtype == modfac.dummy_type:
        typekey = modfac.dummy_gcm_type
    if gcmtype == modfac.ncbased_type:
        typekey = modfac.ncfile_gcm_type
    model = modfac.create_model(typekey, inputdir, workdir,
                                nprocs=gcm_num_procs,
                                redirect=gcm_redirect,
                                channel_type=channel_type,
                                restart=restart,
                                restart_steps=gcm_steps)
    model.initialize_code()
    model.exp_name = gcm_exp_name
    model.num_steps = gcm_steps
    model.step = 0
    model.commit_parameters()
    model.commit_grid()

    log.info("gcm_init called with couple_surface = " + str(couple_surface))
    model.set_vdf_in_sp_mask(not couple_surface)

    return model


# Creates and initialized a LES model
def les_init(lestype, inputdir, workdir, starttime, index):
    typekey = lestype
    if lestype == modfac.dummy_type:
        typekey = modfac.dummy_les_type
    if lestype == modfac.ncbased_type:
        typekey = modfac.ncfile_les_type

    # optionally, schedule restart files to be written at the end of the run
    # currently, one restart is explicitly requested at the end of the run, in filanlize()
    # trestart = gcm_steps * (900 | units.s) # TODO: gcm dt currently hardcoded
    # if not restart:
    #    trestart += les_spinup | units.s             
    #    # add spinup time to trestart

    trestart = 0 | units.s # don't write periodic restarts
    
    model = modfac.create_model(typekey, inputdir, workdir,
                                nprocs=les_num_procs,
                                redirect=les_redirect,
                                channel_type=channel_type,
                                restart=restart,
                                trestart=trestart,
                                starttime=starttime,
                                index=index,
                                qt_forcing=qt_forcing)
#    model.initialize_code()
    model.commit_parameters()
    model.commit_grid()
    return model


# Starts a number of threads for running les models
def start_worker_threads(num_threads):
    worker_threads = []
    work_queue = Queue()
    for i in range(num_threads):
        t = threading.Thread(target=worker, args=(work_queue, i), name="worker " + str(i))
        worker_threads.append(t)
        t.start()
    return work_queue, worker_threads


# a worker thread - using a work queue
def worker(work_queue, i):
    global errorFlag
    while True:
        les, model_time, offset = work_queue.get()
        if les is None:
            log.info("Worker thread %d exiting" % i)
            work_queue.put((None, None))  # put the special quit work back into the queue
            return  # stop this thread
        log.info("Worker thread %d evolves les at index %d to time %s" % (i, les.grid_index, model_time))
        step_les(les, model_time, offset)
        work_queue.task_done()
        log.info("Worker thread %d is done." % i)


# Signals and waits for all worker threads to stop.
def stop_worker_threads(work_queue, worker_threads):
    log.info("Signalling worker threads to quit...")
    work_queue.put((None, None))  # special work - signals the worker_threads to quit.
    log.info("Waiting for worker threads to quit...")
    for w in worker_threads:  # wait for the worker threads to quit
        w.join()


# TODO: with the queue, could let dales instances run immediately when their forcings are set
# instead of setting all forcings then running all
def step_les_models(model_time, work_queue, offset=les_spinup):
    global errorFlag
    les_wall_times = []
    les_profiles = {}
    if not any(les_models):
        return les_wall_times
    if les_queue_threads >= len(les_models):  # Step all dales models in parallel
        if async_evolve:  # evolve all dales models with asynchronous Amuse calls
            reqs = []
            profile_reqs={}
            pool = AsyncRequestsPool()
	    for les in les_models:
               req=les.evolve_model.asynchronous(model_time + (offset | units.s), exactEnd=True)
	       reqs.append(req)
               pool.add_request(req)
               profiles=spcpl.get_les_profiles(les, True)
               profile_reqs[les] = profiles
	       for r in profiles.values():
                   pool.add_request(r)
            # now while the dales threads are working, sync the netcdf to disk
            spio.sync_root()
	    # wait for all threads
            #return pool, requests
            pool.waitall()
            for les, prof in profile_reqs.items():
                les_profiles[les] = {v : r.result() for v,r in prof.items()}
            try:
                les_wall_times = [r.result().value_in(units.s) for r in reqs]
                log.info("async step_les_models() done. Elapsed times:" + str(['%5.1f' % t for t in les_wall_times]))
            except Exception as e:
                log.error("Exception caught while gathering results: %s" % e.message)
    else:  # sequential version
        for les in les_models:
            walltime=step_les(les, model_time, offset)
            les_wal_times.append(walltime)
            profiles=spcpl.get_les_profiles(les, False)
            les_profiles[les] = profiles
        log.info("sequential step_les_models() done. Elapsed times:" + str(['%5.1f' % t for t in les_wall_times]))

    return les_wall_times, les_profiles


# step a dales instance to a given Time
def step_les(les, stoptime, offset=0):
    start = time.time()
    # small les time steps for cloud field gathering
    step_dt = les_dt | units.s
    epsilon = 1 | units.s  # tolerance for fp comparison
    if not les_dt > 0:
        # simply step until caught up
        start = time.time()
        les.evolve_model(stoptime + (offset | units.s), exactEnd=1)
        walltime = start - time.time() 
    else:
        # fixed-length stepping intervals to save statistics during the les run
        start = time.time()
        t = les.get_model_time()
        log.info("Les at point %d starts at %.0f, should reach %.0f with offset %d" % (
            les.grid_index, t.value_in(units.s), stoptime.value_in(units.s), offset))
        while t < stoptime - epsilon + (offset | units.s):
            t += step_dt
            les.evolve_model(t, exactEnd=1)
        walltime = start - time.time()
    return walltime
