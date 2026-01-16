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


# FUNCTION parse_discovery_positions()
# Parses discovery output file and extracts CAAS position numbers

def parse_discovery_positions(discovery_file, genename):
    """
    Parse discovery output file and extract positions that have CAAS.
    
    Args:
        discovery_file: Path to discovery output file
        genename: Gene name to match (e.g., "BRCA1")
    
    Returns:
        set of position numbers (as strings) that had CAAS in discovery
    """
    if not exists(discovery_file):
        print(f"Warning: Discovery file {discovery_file} not found. Processing all positions.")
        return None
    
    caas_positions = set()
    
    try:
        with open(discovery_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("Gene"):  # Skip header
                    continue
                
                try:
                    # Discovery output format: Gene\tTrait\tPosition\t...
                    fields = line.split('\t')
                    if len(fields) >= 3:
                        gene = fields[0]
                        position = fields[2]
                        
                        # Match gene name and extract position number
                        if gene == genename:
                            caas_positions.add(position)
                except:
                    continue
        
        return caas_positions if len(caas_positions) > 0 else None
    
    except Exception as e:
        print(f"Warning: Error parsing discovery file: {e}. Processing all positions.")
        return None


# FUNCTION fetch_caas():
# fetches caas per each thing

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

    # Filter for the number of gaps and missing species

    if len(valid_traits) > 0:

        for trait in valid_traits[:]:  # Iterate over a copy to allow removal

            ### GAPS filtering.

            if maxgaps_fg != "NO" and processed_position.trait2gaps_fg[trait] > int(maxgaps_fg):
                if trait in valid_traits:
                    valid_traits.remove(trait)
                    continue

            # FIX: Changed maxgaps_fg to maxgaps_bg
            if maxgaps_bg != "NO" and processed_position.trait2gaps_bg[trait] > int(maxgaps_bg):
                if trait in valid_traits:
                    valid_traits.remove(trait)
                    continue

            if maxgaps_all != "NO" and processed_position.trait2gaps_fg[trait] + processed_position.trait2gaps_bg[trait] > int(maxgaps_all):
                if trait in valid_traits:
                    valid_traits.remove(trait)
                    continue


            ### Missing filtering.

            if maxmiss_fg != "NO" and processed_position.trait2miss_fg[trait] > int(maxmiss_fg):
                if trait in valid_traits:
                    valid_traits.remove(trait)
                    continue

            # FIX: Changed maxmiss_fg to maxmiss_bg
            if maxmiss_bg != "NO" and processed_position.trait2miss_bg[trait] > int(maxmiss_bg):
                if trait in valid_traits:
                    valid_traits.remove(trait)
                    continue

            if maxmiss_all != "NO" and processed_position.trait2miss_fg[trait] + processed_position.trait2miss_bg[trait] > int(maxmiss_all):
                if trait in valid_traits:
                    valid_traits.remove(trait)
                    continue
            
            ### Pair-aware filtering (only in paired mode with miss_pair enabled)
            
            if multiconfig.paired_mode and miss_pair and trait in valid_traits:
                # Check if missing thresholds are equal (either both explicitly set, or both using overall)
                miss_thresholds_equal = False
                if maxmiss_fg != "NO" and maxmiss_bg != "NO" and maxmiss_fg == maxmiss_bg:
                    miss_thresholds_equal = True
                elif maxmiss_fg == "NO" and maxmiss_bg == "NO" and maxmiss_all != "NO":
                    # Both fg and bg use the overall threshold - treat as equal
                    miss_thresholds_equal = True
                
                if miss_thresholds_equal:
                    miss_pairs_fg = set(processed_position.trait2miss_pairs_fg.get(trait, []))
                    miss_pairs_bg = set(processed_position.trait2miss_pairs_bg.get(trait, []))
                    
                    # Check if missing pairs match
                    if miss_pairs_fg != miss_pairs_bg:
                        valid_traits.remove(trait)
                        continue
                
                # Check if gap thresholds are equal (either both explicitly set, or both using overall)
                gap_thresholds_equal = False
                if maxgaps_fg != "NO" and maxgaps_bg != "NO" and maxgaps_fg == maxgaps_bg:
                    gap_thresholds_equal = True
                elif maxgaps_fg == "NO" and maxgaps_bg == "NO" and maxgaps_all != "NO":
                    # Both fg and bg use the overall threshold - treat as equal
                    gap_thresholds_equal = True
                
                if gap_thresholds_equal:
                    gap_pairs_fg = set(processed_position.trait2gap_pairs_fg.get(trait, []))
                    gap_pairs_bg = set(processed_position.trait2gap_pairs_bg.get(trait, []))
                    
                    # Check if gapped pairs match
                    if gap_pairs_fg != gap_pairs_bg:
                        if trait in valid_traits:
                            valid_traits.remove(trait)
                        continue

     
        output_traits = []

        # Filter for pattern
        for x in valid_traits:
            # Get foreground and background species lists (already filtered for ungapped)
            fg_species = processed_position.trait2ungapped_fg[x][:]
            fg_species.sort()
            bg_species = processed_position.trait2ungapped_bg[x][:]
            bg_species.sort()

            # Extract amino acid for each species to create expanded pattern
            aa_tag_fg = "".join([processed_position.d[sp].split("@")[0] for sp in fg_species])
            aa_tag_bg = "".join([processed_position.d[sp].split("@")[0] for sp in bg_species])

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag, multiconfig, processed_position.d, max_conserved, x)

            if check.caas == True and check.pattern in admitted_patterns:
                output_traits.append(x)

    else:
        output_traits = []

    # Return the line
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))

    traitline = ",".join(output_traits)
    empval = str(int(count)/cycles)

    outline = "\t".join([position_name, count, str(cycles), empval])

    return outline
        

# FUNCTION boot_on_single_alignment()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file, miss_pair=False, max_conserved=0, discovery_file=None):

    the_genename = sliced_object.genename
    print("caastools found", resampled_traits.cycles, "resamplings")
    
    # OPTIMIZATION: Filter to only test positions found in discovery
    positions_list = list(sliced_object.d)
    total_positions = len(positions_list)
    
    if discovery_file:
        print(f"\nFiltering positions based on discovery results from: {discovery_file}")
        caas_positions = parse_discovery_positions(discovery_file, the_genename)
        
        if caas_positions:
            # Filter positions: only keep those that match discovery CAAS positions
            filtered_positions = []
            for pos_dict in positions_list:
                # Extract position number from the position dictionary
                # Format is "AA@position_number" in values
                for species, aa_info in pos_dict.items():
                    pos_num = aa_info.split("@")[1]
                    if pos_num in caas_positions:
                        filtered_positions.append(pos_dict)
                        break
            
            positions_list = filtered_positions
            speedup = total_positions / max(1, len(positions_list))
            print(f"✓ Optimization: Testing only {len(positions_list)} CAAS positions (from {total_positions} total)")
            print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
            print(f"✓ Tests reduced: {total_positions * resampled_traits.cycles:,} → {len(positions_list) * resampled_traits.cycles:,}\n")
        else:
            print("No CAAS positions found in discovery file. Processing all positions.\n")

    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = resampled_traits, species_in_alignment = sliced_object.species), positions_list)

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
    
    print(f"Results written to {output_file}")

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
