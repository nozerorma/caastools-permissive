#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/

__version__ = "2.0.0-paired"

'''
A Convergent Amino Acid Substitution identification 
and analysis toolbox

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).

Pair-aware implementation: Miguel Ramon (miguel.ramon@upf.edu)

MODULE NAME: boot.py
DESCRIPTION: bootstrap function
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: ct

'''


from modules.init_bootstrap import *
from modules.disco import process_position
from modules.caas_id import iscaas
from modules.alimport import *

from os.path import exists
import functools

# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True

    all_g = gfg + gbg
    
    if max_all != "NO" and all_g > int(max_all):
        out = False

    elif max_fg != "NO" and gfg > int(max_fg):
        out = False

    elif max_bg != "NO" and gbg > int(max_bg):
        out = False

    return out


# FUNCTION filter_for_missing()
# filters a trait for its gaps

def filter_for_missings(max_m_bg, max_m_fg, max_m_all, mfg, mbg):

    out = True

    all_m = mfg + mbg
    
    if max_m_all != "NO" and all_m > int(max_m_all):
        out = False

    elif max_m_fg != "NO" and mfg > int(max_m_fg):
        out = False

    elif max_m_bg != "NO" and mbg > int(max_m_bg):
        out = False

    return out



def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, cycles, multiconfig, miss_pair=False, max_conserved=0, admitted_patterns=["1","2","3"]):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    def filter_trait(trait, the_processed_position, max_bg_gaps, max_fg_gaps, max_all_gaps, max_bg_miss, max_fg_miss, max_all_miss, mc, mp, mc_conserved):

        ### GAPS filtering.

        if max_fg_gaps != "NO" and the_processed_position.trait2gaps_fg[trait] > int(max_fg_gaps):
            return False
        # FIX: Changed maxgaps_fg to max_bg_gaps
        if max_bg_gaps != "NO" and the_processed_position.trait2gaps_bg[trait] > int(max_bg_gaps):
            return False
        if max_all_gaps != "NO" and the_processed_position.trait2gaps_fg[trait] + the_processed_position.trait2gaps_bg[trait] > int(max_all_gaps):
           return False

        ### Missings filtering.
        
        if max_fg_miss != "NO" and the_processed_position.trait2miss_fg[trait] > int(max_fg_miss):
            return False
        # FIX: Changed maxmiss_fg to max_bg_miss
        if max_bg_miss != "NO" and the_processed_position.trait2miss_bg[trait] > int(max_bg_miss):
            return False
        if max_all_miss != "NO" and the_processed_position.trait2miss_fg[trait] + the_processed_position.trait2miss_bg[trait] > int(max_all_miss):
            return False
        
        ### Pair-aware filtering (only in paired mode with miss_pair enabled)
        
        if mc.paired_mode and mp:
            # Check if missing thresholds are equal
            if max_fg_miss != "NO" and max_bg_miss != "NO" and max_fg_miss == max_bg_miss:
                miss_pairs_fg = set(the_processed_position.trait2miss_pairs_fg.get(trait, []))
                miss_pairs_bg = set(the_processed_position.trait2miss_pairs_bg.get(trait, []))
                
                # Check if missing pairs match
                if miss_pairs_fg != miss_pairs_bg:
                    return False
            
            # Check if gap thresholds are equal
            if max_fg_gaps != "NO" and max_bg_gaps != "NO" and max_fg_gaps == max_bg_gaps:
                gap_pairs_fg = set(the_processed_position.trait2gap_pairs_fg.get(trait, []))
                gap_pairs_bg = set(the_processed_position.trait2gap_pairs_bg.get(trait, []))
                
                # Check if gapped pairs match
                if gap_pairs_fg != gap_pairs_bg:
                    return False
        
        # Pattern check (always performed after filtering)
        aa_tag_fg_list = the_processed_position.trait2aas_fg[trait]
        aa_tag_fg_list.sort()
        aa_tag_fg = "".join(aa_tag_fg_list)

        aa_tag_bg_list = the_processed_position.trait2aas_bg[trait]
        aa_tag_bg_list.sort()
        aa_tag_bg = "".join(aa_tag_bg_list)

        tag = "/".join([aa_tag_fg, aa_tag_bg])

        check = iscaas(tag, mc, the_processed_position.d, mc_conserved)

        if check.caas == True and check.pattern in admitted_patterns:
            return True
        
        return False
    
    output_traits = filter(functools.partial(
        filter_trait,
        the_processed_position = processed_position,
        max_bg_gaps = maxgaps_bg,
        max_fg_gaps = maxgaps_fg,
        max_all_gaps = maxgaps_all,

        max_bg_miss = maxmiss_bg,
        max_fg_miss = maxmiss_fg,
        max_all_miss = maxmiss_all,
        
        mc = multiconfig,
        mp = miss_pair,
        mc_conserved = max_conserved,
    ),
                    valid_traits)



    output_traits = list(output_traits)

    # Return the line
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))

    # traitline = ",".join(output_traits) NOT USED ANYMORE
    empval = str(int(count)/cycles)

    outline = "\t".join([position_name, count, str(cycles), empval])

    return outline
        

# FUNCTION boot_on_single_alignment()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file, miss_pair=False, max_conserved=0):


    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = resampled_traits, species_in_alignment = sliced_object.species), sliced_object.d)
    the_genename = sliced_object.genename
    print("caastools found", resampled_traits.cycles, "resamplings")

    # Step 4: extract the raw caas

    output_lines = map(
        functools.partial(
            caasboot,
            list_of_traits = resampled_traits.alltraits,
            genename = the_genename,
            maxgaps_fg = max_fg_gaps,
            maxgaps_bg = max_bg_gaps,
            maxgaps_all = max_overall_gaps,

            maxmiss_fg = max_fg_miss,
            maxmiss_bg = max_bg_miss,
            maxmiss_all = max_overall_miss,

            multiconfig = resampled_traits,
            miss_pair = miss_pair,
            max_conserved = max_conserved,

            admitted_patterns = the_admitted_patterns,
            cycles = resampled_traits.cycles) ,processed_positions
    )

    output_lines = list(output_lines)

    ooout = open(output_file, "w")

    for line in output_lines:
        print(line, file=ooout)
    
    ooout.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(bootstrap_result):
    with open(bootstrap_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d
