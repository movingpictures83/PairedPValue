
import pandas as pd
import numpy as np
from scipy import stats
import os

def get_diff(row):
    return row["PTR_y"] - row["PTR_x"]

import PyIO
import PyPluMA

class PairedPValuePlugin:
 def input(self, inputfile):
   self.parameters = PyIO.readParameters(inputfile)

 def run(self):
     pass

 def output(self, outputfile):
   level = self.parameters["level"]#"species"
   out_file = outputfile#"Significance.txt"

   ptr_df = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["ptr"])#pd.read_csv("all_ptr.csv")
   mylist = PyIO.readSequential(PyPluMA.prefix()+"/"+self.parameters["columns"])
   ptr_df = ptr_df[mylist]
   #ptr_df = ptr_df[["Sample_ID", "coverage", "genome", "PTR", "Day_of_Life",
   #              "Individual", "Antibiotic_Treatment", "species", "genus", "family", "class", "order", "phylum"]]
   treatments = list(ptr_df["Antibiotic_Treatment"].unique())

   treatments.remove(np.nan)

   treatment_names = [x.split("-Before")[0] for x in treatments]
   treatment_names = [x.split("-After")[0] for x in treatments]

   treatment_names = list(set(treatment_names))

   # For each treatment, identify paired p-value of each strain weather it has statistically significant from before and after

   #if os.path.exists(out_file):
   #    os.remove(out_file)

   with open(out_file, "w") as out:
    out.write("species\tchange\tp_value\tSE\tTreatment\n")

    for treatment in treatment_names:
        # Get the list of before
        before_df = ptr_df[ptr_df["Antibiotic_Treatment"]==treatment+"-Before"].reset_index(drop=True)
        after_df = ptr_df[ptr_df["Antibiotic_Treatment"]==treatment+"-After"].reset_index(drop=True)

        # Merge before and after
        merged_df = before_df.merge(after_df, how="outer", on=["genome", "Individual"])
        merged_df = merged_df.dropna()

        unique_species = list(merged_df[level + "_x"].unique())
        if len(unique_species)>0:
            #print(".. Analyzing treatment {}".format(treatment))
            #out.write(treatment+"\n")

            for species in unique_species:
                df_species = merged_df[merged_df["species_x"]==species]
                df_species["diff"] = df_species.apply(get_diff, axis=1)
                p_value = stats.ttest_rel(df_species["PTR_x"], df_species["PTR_y"]).pvalue
                if species=="Serratia sp. FDAARGOS_506":
                    pass
                    pass
                if not np.isnan(p_value):
                    difference = df_species["diff"].mean()
                    SE = df_species["diff"].sem()
                    out.write("{}\t{}\t{}\t{}\t{}\n".format(species, difference, p_value, SE, treatment))

                    #print("{}: p_value: {}, difference: {}".format(species, p_value, difference))




    pass
    pass
