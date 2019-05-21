
import tables as tb
import reco_functions as rf
from components import part_wfs_from_files
from components import true_photoelect
from components import wvf_passing_filter

from invisible_cities.reco                import tbl_functions as tbl
from invisible_cities.dataflow            import dataflow      as fl
from invisible_cities.dataflow.dataflow   import push
from invisible_cities.dataflow.dataflow   import pipe
from invisible_cities.io.run_and_event_io import run_and_event_writer

from invisible_cities.cities.components import city
from invisible_cities.cities.components import print_every
from invisible_cities.cities.components import compute_xy_position


@city
def city_r_maps(files_in, file_out, compression, event_range, print_mod):

    #sens_pos       = rf.sensor_position    (h5in)
    #sens_pos_cyl   = rf.sensor_position_cyl(h5in)

    true_phot = fl.map(true_photoelect(compton=False),
                       args="particles",
                       out=("true_pos1", "true_pos2"))

    wvf_pass_filter = fl.map(wvf_passing_filter(nsteps, th_start),
                   args=("wfs", "true_pos1", "true_pos2", "sens_pos", "sens_pos_cyl"),
                   out=("region1_info", "region2_info"))


    event_count_in        = fl.spy_count()
    event_count_out       = fl.spy_count()

    with tb.open_file(file_out, "w", filters = tbl.filters(compression)) as h5out:

        write_event_info      = fl.sink(run_and_event_writer(h5out                ), args=("run_number", "event_number", "timestamp"))

        return push(source = part_wfs_from_files(files_in), ,
                    pipe   = pipe(
                        fl.slice(*event_range, close_all=True),
                        print_every(print_mod)                ,
                        event_count_in       .spy             ,
                        true_phot                             ,
                        wvf_pass_filter                       ,
                        event_count_out      .spy             ,
                        fl.fork(write_event_info              )),
                    result = dict(events_in  = event_count_in .future,
                                  events_out = event_count_out.future))
