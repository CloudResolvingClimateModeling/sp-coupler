# from amuse.rfi import channel
import logging
import numpy

# Logger
log = logging.getLogger(__name__)


# This method broadcasts the process types to all MPI tasks,
# allowing them to set up model communicators
def send_model_colors(gcm_procs, les_procs, num_les):
    from mpi4py import MPI
    worldcomm = MPI.COMM_WORLD
    procs = worldcomm.Get_size()
    if 1 + gcm_procs + num_les * les_procs > procs:
        raise Exception("Too few processes launched for requested model configuration")
    if 1 + gcm_procs + num_les * les_procs < procs:
        log.warning("Too many processes launched for requested model configuration")
    # We assume the process has been launched as mpiexec -n1 spmaster.py -nx gcm_worker -ny les_worker
    # colors = numpy.repeat(-1,procs)
    colors = numpy.zeros(procs, dtype=numpy.int32)
    colors[0] = 1
    colors[1:gcm_procs + 1] = 2
    for i in range(num_les):
        offset = 1 + gcm_procs + i * les_procs
        colors[offset:offset + les_procs] = 3 + i
    log.info("Master process scattering color array " + str(colors))

    rec = numpy.zeros(procs, dtype=numpy.int32)
    worldcomm.Scatter([colors, 1, MPI.INT32_T], [rec, 1, MPI.INT32_T], root=0)
    log.info("Master process received " + str(rec))

    worldcomm.Split(colors[0], 0)


    # moved this into amuse/src/amuse/rfi/channel.py
# class MpiChannelNoSpawn(channel.MpiChannel):

#     process_counter = 1
#     code_communicators = [MPI.COMM_SELF]

#     def start(self):
#         if not self.debugger_method is None:
#             command,arguments = self.debugger_method(self.full_name_of_the_worker,
#                                                      self,
#                                                      interpreter_executable = self.interpreter_executable)
#         else:
#             if not self.can_redirect_output or
#                   (self.redirect_stdout_file == 'none' and self.redirect_stderr_file == 'none'):

#                 if self.interpreter_executable is None:
#                     command = self.full_name_of_the_worker
#                     arguments = None
#                 else:
#                     command = self.interpreter_executable
#                     arguments = [self.full_name_of_the_worker]
#             else:
#                 command,arguments = self.REDIRECT(self.full_name_of_the_worker,
#                                                   self.redirect_stdout_file,
#                                                   self.redirect_stderr_file,
#                                                   command = self.python_exe_for_redirection,
#                                                   interpreter_executable = self.interpreter_executable)

#         log.info("MpiChannelNoSpawn.start")

#         proc_offset = MpiChannelNoSpawn.process_counter
#         worker_processes = range(proc_offset, proc_offset + self.number_of_workers)

#         #group = MPI.COMM_WORLD.Get_group().Incl(worker_processes)
#         #peer_comm = MpiChannelNoSpawn.code_communicators[process_counter]
#        # self.intercomm = MPI.COMM_SELF.Create_intercomm(0,peer_comm,0)
#         # MpiChannelNoSpawn.process_counter += self.number_of_workers
