import tables as tb
import numpy  as np

from antea.io.mc_io                   import read_mcsns_response
from invisible_cities.io  .mcinfo_io  import read_mcinfo
from invisible_cities.core.exceptions import NoHits


def get_true_pos(part):
    hit_positions = [h.pos for h in part.hits]
    energies      = [h.E   for h in part.hits]
    energy        = sum(energies)
    if energy: return np.average(hit_positions, axis=0, weights=energies)
    else: raise NoHits


def true_photoelect(h5in, true_file, evt, compton=False):

    this_event_dict = read_mcinfo        (     h5in, (evt, evt+1))
    this_event_wvf  = read_mcsns_response(true_file, (evt, evt+1))
    part_dict       = list(this_event_dict.values())[0]

    ave_true1 = []
    ave_true2 = []

    for indx, part in part_dict.items():
        if part.name == 'e-' :
            mother = part_dict[part.mother_indx]
            if part.initial_volume == 'ACTIVE' and part.final_volume == 'ACTIVE':
                if mother.primary and np.isclose(mother.E*1000., 510.999, atol=1.e-3):
                    if compton==True: pass
                    else:
                        if np.isclose(sum(h.E for h in part.hits), 0.476443, atol=1.e-6): pass
                        else: continue

                    if mother.p[1] > 0.: ave_true1 = get_true_pos(part)
                    else:                ave_true2 = get_true_pos(part)
    return ave_true1, ave_true2


def find_closest_sipm(x, y, z, sens_pos, sns_over_thr, charges_over_thr):
    """Returns the sensor ID of the closest sipm to the given true event
    """
    min_dist = 1.e9
    min_sns  = 0
    for sns_id, charge in zip(sns_over_thr, charges_over_thr):
        pos  = sens_pos[sns_id]
        dist = np.linalg.norm(np.subtract((x, y, z), pos))
        if dist < min_dist:
            min_dist = dist
            min_sns  = sns_id
    return min_sns
