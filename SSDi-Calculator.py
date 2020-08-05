import argparse
import logging
import os
import random
import numpy as np
from scipy import stats
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
      SSDi-Calculator.py: Calculate sexual size dimorphism index and associated statistics. 
      ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--input",
                            required=True,
                            help="The full path to a text file containing "
                            "the input data.")
    
    parser.add_argument("-f", "--fileformat",
                            required=True,
                            choices=["tab","csv"],
                            help="The format of the input text file, either "
                            "tab-delimited (tab) or comma-separated values (csv).")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="The full path to an existing directory "
                            "to write output all files.")
    
    return parser.parse_args()

def check_inputs(input, outdir):
    """
    Function to check that input file paths and output 
    directory are valid. If not, raises IOError.

    Args:
    input - Full path to the input data file.
    outdir - Full path to an existing output directory.
    """

    # make sure transcipt file path is valid
    if not os.path.isfile(input):
        raise IOError("\n\n{} does not appear to be a valid file.".format(input))
    
    # check output directory path is valid
    if not os.path.isdir(outdir):
        raise IOError("\n\n{} is not a valid directory.".format(outdir))

def setup_logging():
    """
    Set up logging to file and on-screen.
    """
    # ensure logging file does not exist, if so remove
    if os.path.exists("SSDi-Calculator-Run.log"):
        os.remove("SSDi-Calculator-Run.log")
        
    # set up logging to file
    logging.basicConfig(filename="SSDi-Calculator-Run.log",
                            format="%(levelname)s: %(asctime)s: %(message)s",
                            datefmt='%d-%b-%y %H:%M:%S',
                            level=logging.DEBUG)
    
    # set up logging to console 
    console = logging.StreamHandler() 
    console.setLevel(logging.INFO) 
    formatter = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s",
                                      datefmt='%H:%M:%S')
    console.setFormatter(formatter) 
    logging.getLogger().addHandler(console) 

def input_to_dict(input, fileformat):
    """
    Build a dictionary data structure from input file,
    where each key is a species and the value is a sub-dictionary 
    of sex-specific sizes. The sub-dictionary contains keys "M" and "F", 
    with values of lists containing floats of size data. Returns the 
    dictionary, labeled as species_dict.

    Args:
    input - Full path to the input data file.
    fileformat - Specifies tab-delimited or csv format for input file. 
    """    
    print("\n\n\n")
    logging.info("Reading input data file: {}".format(input.split('/')[-1]))
    
    # initiate empty dictionary
    species_dict = {}

    # open transcript file
    with open(input, 'r') as fh:
        # skip header
        next(fh)
        # iterate over lines in file
        for line in fh:
            # if line is not blank
            if line.strip():
                # split line by tab or column, depending on format specified
                if fileformat == 'tab':
                    cols = line.strip().split('\t')
                elif fileformat == 'csv':
                    cols = line.strip().split(',')

                # quick check that expected number of columns is present
                if len(cols) < 3:
                    logging.debug("input_to_dict: Bad line - {}".format(line))
                # if 3 or more cols:
                else:
                    # label split elements for simplicity
                    species, sex, size = cols[0], cols[1].upper(), float(cols[2])

                    # check if key is in dictionary
                    if species not in species_dict:
                        # if not, add species key to dict, with new sex-specific dict
                        # with keys = M, F, and vals as empty lists
                        species_dict[species] = {"M":[], "F":[]}
                        # add data
                        if sex == "M":
                            species_dict[species]["M"].append(size)
                        elif sex == "F":
                            species_dict[species]["F"].append(size)
                        else:
                            logging.debug("input_to_dict: Skipping entry - {}".format(species, sex, size))
                    else:
                        # if key already in dictionary, add data
                        if sex == "M":
                            species_dict[species]["M"].append(size)
                        elif sex == "F":
                            species_dict[species]["F"].append(size)
                        else:
                            logging.debug("input_to_dict: Skipping entry - {}, {}, {}".format(species, sex, size))
                            
    logging.info("Found data for {} species.\n\n".format(len(species_dict)))
    
    # send back data dictionary
    return species_dict

def quick_counts(species_dict):
    """
    Log information about sample sizes for each species.

    Args:
    species_dict - Dictionary data structure from input_to_dict().
    """

    # initiate three different zero counts
    onecount, twocount, threecount = int(0), int(0), int(0)

    # iterate over dictionary keys
    for k, v in sorted(species_dict.items()):
        #print("{}:\n\tM:{}\n\tF:{}\n".format(k, v["M"], v["F"]))
        # add to relevant counters
        if len(v["M"]) > 1 and len(v["F"]) > 1:
            onecount += 1
        if len(v["M"]) > 2 and len(v["F"]) > 2:
            twocount += 1
        if len(v["M"]) > 3 and len(v["F"]) > 3:
            threecount += 1
            
    # log information
    logging.debug("Species with > 1 M and > 1 F: {}".format(onecount))
    logging.debug("Species with > 2 M and > 2 F: {}".format(twocount))
    logging.debug("Species with > 3 M and > 3 F: {}\n\n".format(threecount))
    

def run_ssdi_calculations(species_dict):
    """
    Runs all tests using data stored in species_dict 
    data structure.

    Args:
    species_dict - Dictionary data structure from input_to_dict().
    """

    # initiate empty dictionary to store all results
    calc_dict = {}
    
    # iterate over dictionary keys
    for k, v in sorted(species_dict.items()):
        # log species and counts per sex
        logging.info("Species: {}".format(k))
        logging.info("Males: {}".format(len(v["M"])))
        logging.info("Females: {}".format(len(v["F"])))

        # series of rules to handle edge cases properly
        if len(v["F"]) == 1 and len(v["M"]) == 1:
            # run standard ssdi calculation using each point estimate
            ssdi = ssdi_single(v["F"][0], v["M"][0])
            logging.info("Standard SSDi: {}\n".format(ssdi))
            # no pairwise comparisons are possible, skip remaining tests
            avg_ssdi, p_pair, low, high, p_perm, diff, avgf, avgm = "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            
        elif len(v["F"]) > 1 and len(v["M"]) == 1:
            # run standard ssdi calculation using average from females and point estimate for males
            ssdi = ssdi_single(round(np.mean(v["F"]), 1), v["M"][0])
            logging.info("Standard SSDi: {}.".format(ssdi))
            # perform pairwise calculations and corresponding t-test
            avg_ssdi, p_pair = ssdi_pairwise(v["F"], v["M"], ttest=True)
            logging.info("Pairwise Analyses: Average pairwise SSDi: {}.".format(avg_ssdi))
            logging.debug("Pairwise Analyses: One-sample t-test P-value: {}.".format(p_pair))
            # perform permutation test
            low, high, p_perm = run_permutations(v["F"], v["M"], avg_ssdi)
            diff = abs(ssdi - avg_ssdi)
            avgf, avgm = round(np.mean(v["F"]), 1), "NA"

        elif len(v["F"]) == 1 and len(v["M"]) > 1:
            # run standard ssdi calculation using point estimate for females and average from males
            ssdi = ssdi_single(v["F"][0], round(np.mean(v["M"]), 1))
            logging.info("Standard SSDi: {}.".format(ssdi))
            # perform pairwise calculations and corresponding t-test
            avg_ssdi, p_pair = ssdi_pairwise(v["F"], v["M"], ttest=True)
            logging.info("Pairwise Analyses: Average pairwise SSDi: {}.".format(avg_ssdi))
            logging.debug("Pairwise Analyses: One-sample t-test P-value: {}.".format(p_pair))
            # perform permutation test
            low, high, p_perm = run_permutations(v["F"], v["M"], avg_ssdi)
            diff = abs(ssdi - avg_ssdi)
            avgf, avgm = "NA", round(np.mean(v["M"]), 1)

        elif len(v["F"]) > 1 and len(v["M"]) > 1:
            # run standard ssdi calculation using averages per sex
            ssdi = ssdi_single(round(np.mean(v["F"]), 1), round(np.mean(v["M"]), 1))
            logging.info("Standard SSDi: {}.".format(ssdi))
            # perform pairwise calculations and corresponding t-test
            avg_ssdi, p_pair = ssdi_pairwise(v["F"], v["M"], ttest=True)
            logging.info("Pairwise Analyses: Average pairwise SSDi: {}.".format(avg_ssdi))
            logging.debug("Pairwise Analyses: One-sample t-test P-value: {}.".format(p_pair))
            # perform permutation test
            low, high, p_perm = run_permutations(v["F"], v["M"], avg_ssdi)
            diff = abs(ssdi - avg_ssdi)
            avgf, avgm = round(np.mean(v["F"]), 1), round(np.mean(v["M"]), 1)

        else:
            # these are species missing data for one of the sexes
            # we will skip them and not include them in the output file
            ssdi = None
            logging.error("Species {} does not have at least 1 M and 1 F, skipping calculations.\n".format(k))

        # species with data will have ssdi val, use to eliminate bad species
        if ssdi:
            # add all relevant vals to new dictionary structure
            calc_dict[k] = {"males":len(v["M"]), "females":len(v["F"]),
                                "ssdi":ssdi, "avg_ssdi":avg_ssdi, "p_pair":p_pair,
                                "low":low, "high":high, "p_perm":p_perm, "diff":diff,
                                "avgf":avgf, "avgm":avgm}
    return calc_dict

    
def ssdi_single(f, m):
    """
    Function to calculate SSDi from Lovich & Gibbons 1992, where
    SSDi = [(larger sex / smaller sex) - 1], arbitrarily set negative 
    if males are larger and positive if females are larger. Returns 
    the SSDi val.

    Args:
    f - Single value of female size (float)
    m - Single value of male size (float)
    """
    if m > f:
        ssdi = round((-1*((m / f) - 1)), 3)
    elif f > m:
        ssdi = round(((f / m) - 1), 3)
    elif f == m:
        ssdi = float(0)
        
    return ssdi
    

def ssdi_pairwise(females, males, ttest=False):
    """
    Performs all pairwise SSDi calculations using lists of floats 
    for males and females. Optionally computes 1-sample t-test 
    to compare empirical mean to hypothesized mean of 0, and 
    returns p-val. A p-val < 0.05 allows rejection of null, meaning 
    empirical mean is different from 0. Nonsignificant p-val indicates 
    empirical mean is not different from 0.

    Args:
    females - List of float values for female body sizes.
    males - List of float values for female body sizes.
    ttest - If True, runs 1-sample t-test; Default = False.
    """

    # intiate empty list to store pairwise ssdi vals
    vals = []
    # iterate over females
    for f in females:
        # iterate over males
        for m in males:
            # for each combination, calculate ssdi
            vals.append(ssdi_single(f, m))

    # get mean pairwise ssdi val
    avg_ssdi = round(np.mean(vals), 3)

    # if no t-test, assign NA to p-val
    if ttest is False:
        p = "NA"

    else:
        # perform 1 sample t-test against hypothesized mean of 0
        results = stats.ttest_1samp(vals, 0)

        # get p-val in readable format
        if results[1] <= 0.001:
            p = "<0.001"
        else:
            p = round(results[1], 3)

    return avg_ssdi, p

def run_permutations(females, males, emp_ssdi):
    """
    Performs permutation test with 10,000 bootstraps to 
    generate a null distribution where male size = female size. 
    Shuffles male and female labels and performs all pairwise 
    comparisons, then calculates mean SSDi for the replicate.
    Based on resulting distribution from bootstrapping, calculates 
    a two-tailed p-val based on position of empirical SSDi. 
    Returns critical values from bootstrapping distribution and 
    the p-val.

    Args:
    females - List of float values for female body sizes.
    males - List of float values for female body sizes.
    emp_ssdi - Empirical SSDi value.
    """

    # create initial combined list of males + females
    allsizes = females + males
    #print "combined", allsizes
    # initiate empty list to store permuted ssdi means
    permuted_ssdi_vals = []
    # perform 10,000 bootstraps
    for i in range(0, 10000):
        # randomly shuffle combined list
        randomized = random.sample(allsizes, len(allsizes))
        #print "random", randomized
        # create new lists of females and males that are the same
        # length as the originals
        newf, newm = randomized[:len(females)], randomized[len(females):]
        #print "newf", newf
        #print "newm", newm
        # get mean ssdi
        perm_avg_ssdi, p = ssdi_pairwise(newf, newm)
        permuted_ssdi_vals.append(perm_avg_ssdi)

    # perform significance testing
    permuted_ssdi_vals.sort()
    
    # get percentiles for test
    low, high = round(np.percentile(permuted_ssdi_vals, 2.5), 3), round(np.percentile(permuted_ssdi_vals, 97.5), 3)
    
    # calculate p-value based on position of empirical value in distribution
    ptwotail = round((float(len([x for x in permuted_ssdi_vals if abs(x) > abs(emp_ssdi)])) + 1) / (float(len(permuted_ssdi_vals)) + 1), 5)
    # get p-val in readable format
    if ptwotail <= 0.001:
        p = "<0.001"
    else:
        p = round(ptwotail, 3)

    # write critical vals and empirical val to log
    logging.info("Permutation Test: 2.5 and 97.5 percentile values: {}, {}.".format(low, high))
    logging.info("Permutation Test: Empirical value: {}".format(emp_ssdi))

    # log statement about where empirical value is in distribution
    if emp_ssdi <= low:
        logging.info("Permutation Test: Empirical value lies outside the 2.5 percentile.")

    elif emp_ssdi > low and emp_ssdi < high:
        logging.info("Permutation Test: Empirical value within the 2.5 and 97.5 percentiles.")

    elif emp_ssdi >= high:
        logging.info("Permutation Test: Empirical value lies outside the 97.5 percentile.")

    # log p-val of permutation test
    logging.info("Permutation Test: Permutation test P-value: {}.\n".format(p))

    return low, high, p


def write_output(calc_dict):
    """
    Write all information to tab-delimited and csv output files.

    Args:
    calc_dict - Dictionary data structure resulting from run_ssdi_calculations().
    """
    logging.info("Writing output files.\n\n")
    
    # ensure output files do not exist, if so remove
    if os.path.exists("SSDi-Results.txt"):
        os.remove("SSDi-Results.txt")
    if os.path.exists("SSDi-Results.csv"):
        os.remove("SSDi-Results.csv")

    # write all relevant data to outputs
    with open('SSDi-Results.txt', 'a') as fh1, open('SSDi-Results.csv', 'a') as fh2:
        
        fh1.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format("Species", "Number_Males", "Number_Females",
                                                                                        "Avg_Male", "Avg_Female", "Standard_SSDi",
                                                                                        "Avg_Pairwise_SSDi", "AbsDifference",
                                                                                        "Dimorphism_PValue", "2.5_percentile",
                                                                                        "97.5_percentile"))
        
        fh2.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n'.format("Species", "Number_Males", "Number_Females",
                                                                              "Avg_Male", "Avg_Female","Standard_SSDi",
                                                                              "Avg_Pairwise_SSDi", "AbsDifference",
                                                                              "Dimorphism_PValue", "2.5_percentile",
                                                                              "97.5_percentile"))
        
        for k, v in sorted(calc_dict.items()):
            
            fh1.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(k, calc_dict[k]["males"], calc_dict[k]["females"],
                                                                                            calc_dict[k]["avgm"], calc_dict[k]["avgf"],
                                                                                            calc_dict[k]["ssdi"], calc_dict[k]["avg_ssdi"],
                                                                                            calc_dict[k]["diff"], calc_dict[k]["p_perm"],
                                                                                            calc_dict[k]["low"], calc_dict[k]["high"]))
            
            fh2.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n'.format(k, calc_dict[k]["males"], calc_dict[k]["females"],
                                                                                  calc_dict[k]["avgm"], calc_dict[k]["avgf"],
                                                                                  calc_dict[k]["ssdi"], calc_dict[k]["avg_ssdi"],
                                                                                  calc_dict[k]["diff"], calc_dict[k]["p_perm"],
                                                                                  calc_dict[k]["low"], calc_dict[k]["high"]))
    
    
def main():
    tb = datetime.now()
    args = get_args()
    
    # ensure CL arguments are valid files/paths
    check_inputs(args.input, args.outdir)

    # move to output directory
    os.chdir(args.outdir)

    # configure logging
    setup_logging()
    logging.debug("BEGINNING ANALYSIS")
    logging.debug("Arguments: \n\t\t-i {} \n\t\t-f {} \n\t\t-o {}".format(args.input, args.fileformat, args.outdir))
    
    # get information from transcripts input file
    species_dict = input_to_dict(args.input, args.fileformat)

    # record quick info about sample sizes
    quick_counts(species_dict)

    # create data structure to store all test results
    calc_dict = run_ssdi_calculations(species_dict)

    # write output files
    write_output(calc_dict)
            
    logging.info("Finished! Total elapsed time: {} (H:M:S)\n\n".format(datetime.now() - tb))

    
if __name__ == '__main__':
    main()

