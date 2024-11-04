# This file contains the settings used to measure performance, including
# time spent in certain critical MPI functions.
#
# This file can be `source`d before running `appositionizer` to activate
# caliper measurements.

export CALI_CHANNEL_FLUSH_ON_EXIT=false
export CALI_MPIREPORT_FILENAME=caliper.json
export CALI_MPIREPORT_LOCAL_CONFIG="SELECT sum(sum#time.duration),
                                           inclusive_sum(sum#time.duration)
                                    GROUP BY prop:nested"
export CALI_MPIREPORT_CONFIG='SELECT annotation,
                                     mpi.function,
                                     min(sum#sum#time.duration) as "exclusive_time_rank_min",
                                     max(sum#sum#time.duration) as "exclusive_time_rank_max",
                                     avg(sum#sum#time.duration) as "exclusive_time_rank_avg",
                                     min(inclusive#sum#time.duration) AS "inclusive_time_rank_min",
                                     max(inclusive#sum#time.duration) AS "inclusive_time_rank_max",
                                     avg(inclusive#sum#time.duration) AS "inclusive_time_rank_avg",
                                     percent_total(sum#sum#time.duration) AS "Exclusive time %",
                                     inclusive_percent_total(sum#sum#time.duration) AS "Inclusive time %"
                                     GROUP BY prop:nested FORMAT json'
export CALI_SERVICES_ENABLE=aggregate,event,mpi,mpireport,timestamp
# Everything not blacklisted is profiled. This list was stolen from Caliper...
export CALI_MPI_BLACKLIST="MPI_Comm_rank,MPI_Comm_size,MPI_Wtick,MPI_Wtime"
