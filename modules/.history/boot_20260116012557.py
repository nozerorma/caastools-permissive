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
import numpy as np
from scipy import sparse
import sys
import time
from multiprocessing import Pool, cpu_count

# Global variables for multiprocessing worker processes
_worker_state = {}

def _init_worker(resampled_traits, species_in_alignment, alltraits, genename, 
                 max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, 
                 max_bg_miss, max_overall_miss, miss_pair, max_conserved, 
                 admitted_patterns, cycles):
    """Initialize worker process with shared state"""
    global _worker_state
    _worker_state = {
        'resampled_traits': resampled_traits,
        'species_in_alignment': species_in_alignment,
        'alltraits': alltraits,
        'genename': genename,
        'max_fg_gaps': max_fg_gaps,
        'max_bg_gaps': max_bg_gaps,
        'max_overall_gaps': max_overall_gaps,
        'max_fg_miss': max_fg_miss,
        'max_bg_miss': max_bg_miss,
        'max_overall_miss': max_overall_miss,
        'miss_pair': miss_pair,
        'max_conserved': max_conserved,
        'admitted_patterns': admitted_patterns,
        'cycles': cycles
    }


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
            # Check if missing thresholds are equal (either both explicitly set, or both using overall)
            miss_thresholds_equal = False
            if max_fg_miss != "NO" and max_bg_miss != "NO" and max_fg_miss == max_bg_miss:
                miss_thresholds_equal = True
            elif max_fg_miss == "NO" and max_bg_miss == "NO" and max_all_miss != "NO":
                # Both fg and bg use the overall threshold - treat as equal
                miss_thresholds_equal = True
            
            if miss_thresholds_equal:
                miss_pairs_fg = set(the_processed_position.trait2miss_pairs_fg.get(trait, []))
                miss_pairs_bg = set(the_processed_position.trait2miss_pairs_bg.get(trait, []))
                
                # Check if missing pairs match
                if miss_pairs_fg != miss_pairs_bg:
                    return False
            
            # Check if gap thresholds are equal (either both explicitly set, or both using overall)
            gap_thresholds_equal = False
            if max_fg_gaps != "NO" and max_bg_gaps != "NO" and max_fg_gaps == max_bg_gaps:
                gap_thresholds_equal = True
            elif max_fg_gaps == "NO" and max_bg_gaps == "NO" and max_all_gaps != "NO":
                # Both fg and bg use the overall threshold - treat as equal
                gap_thresholds_equal = True
            
            if gap_thresholds_equal:
                gap_pairs_fg = set(the_processed_position.trait2gap_pairs_fg.get(trait, []))
                gap_pairs_bg = set(the_processed_position.trait2gap_pairs_bg.get(trait, []))
                
                # Check if gapped pairs match
                if gap_pairs_fg != gap_pairs_bg:
                    return False
        
        # Pattern check (always performed after filtering)
        # Get foreground and background species lists (already filtered for ungapped)
        fg_species = the_processed_position.trait2ungapped_fg[trait][:]
        fg_species.sort()
        bg_species = the_processed_position.trait2ungapped_bg[trait][:]
        bg_species.sort()

        # Extract amino acid for each species to create expanded pattern
        aa_tag_fg = "".join([the_processed_position.d[sp].split("@")[0] for sp in fg_species])
        aa_tag_bg = "".join([the_processed_position.d[sp].split("@")[0] for sp in bg_species])

        tag = "/".join([aa_tag_fg, aa_tag_bg])

        check = iscaas(tag, mc, the_processed_position.d, mc_conserved, trait)

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


# FUNCTION _process_single_position_wrapper()
# Module-level wrapper for multiprocessing (nested functions can't be pickled)

def _process_single_position_wrapper(pos_dict, resampled_traits, species_in_alignment, alltraits, genename, 
                                      max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, 
                                      max_overall_miss, miss_pair, max_conserved, admitted_patterns, cycles):
    """Process a single position through the full pipeline (module-level for pickling)"""
    # Step 1: Process position
    processed = process_position(
        pos_dict,
        multiconfig=resampled_traits,
        species_in_alignment=species_in_alignment
    )
    
    # Step 2: Apply bootstrap analysis
    result = caasboot(
        processed,
        list_of_traits=alltraits,
        genename=genename,
        maxgaps_fg=max_fg_gaps,
        maxgaps_bg=max_bg_gaps,
        maxgaps_all=max_overall_gaps,
        maxmiss_fg=max_fg_miss,
        maxmiss_bg=max_bg_miss,
        maxmiss_all=max_overall_miss,
        multiconfig=resampled_traits,
        miss_pair=miss_pair,
        max_conserved=max_conserved,
        admitted_patterns=admitted_patterns,
        cycles=cycles
    )
    
    return result
        

# FUNCTION boot_on_single_alignment()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue
# OPTIMIZED VERSION with vectorization, multiprocessing, progress reporting, and selective position filtering

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file, miss_pair=False, max_conserved=0, discovery_file=None, n_jobs=None):

    the_genename = sliced_object.genename
    print("caastools found", resampled_traits.cycles, "resamplings")
    
    # Ensure sparse matrices are built
    if not resampled_traits.matrix_built:
        print("Building sparse matrices for vectorized processing...")
        resampled_traits.build_sparse_matrices()
    
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
            print("No CAAS positions found in discovery file. Testing all positions.\n")
    
    n_positions = len(positions_list)
    
    # Determine number of parallel workers
    if n_jobs is None:
        n_jobs = min(cpu_count(), n_positions)
    else:
        n_jobs = min(n_jobs, cpu_count(), n_positions)
    
    print(f"Processing {n_positions} positions using {n_jobs} parallel workers...")
    
    # Create partial function with all fixed parameters
    process_func = functools.partial(
        _process_single_position_wrapper,
        resampled_traits=resampled_traits,
        species_in_alignment=sliced_object.species,
        alltraits=resampled_traits.alltraits,
        genename=the_genename,
        max_fg_gaps=max_fg_gaps,
        max_bg_gaps=max_bg_gaps,
        max_overall_gaps=max_overall_gaps,
        max_fg_miss=max_fg_miss,
        max_bg_miss=max_bg_miss,
        max_overall_miss=max_overall_miss,
        miss_pair=miss_pair,
        max_conserved=max_conserved,
        admitted_patterns=the_admitted_patterns,
        cycles=resampled_traits.cycles
    )
    
    # Open output file
    ooout = open(output_file, "w")
    
    # Process positions in parallel with progress tracking
    start_time = time.time()
    positions_processed = 0
    
    print("Progress: ", end="", flush=True)
    
    # Use multiprocessing Pool for parallel execution
    with Pool(processes=n_jobs) as pool:
        # imap returns results as they complete (preserves order)
        for result_line in pool.imap(process_func, positions_list):
            print(result_line, file=ooout)
            positions_processed += 1
            
            # Update progress bar
            progress_pct = int((positions_processed / n_positions) * 100)
            progress_bar_len = int((positions_processed / n_positions) * 50)
            sys.stdout.write(f"\rProgress: [{'=' * progress_bar_len}{' ' * (50 - progress_bar_len)}] {progress_pct}% ({positions_processed}/{n_positions})")
            sys.stdout.flush()
    
    ooout.close()
    
    elapsed_time = time.time() - start_time
    print(f"\nCompleted in {elapsed_time:.2f} seconds ({positions_processed/elapsed_time:.1f} positions/sec)")
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
