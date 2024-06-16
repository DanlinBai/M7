# If one is interested in other reactions (e.g., glycolysis, oxphos, etc.), 
# one needs to search for Reaction_ID using the human_gem file. 
data("human_gem")
data("human_blood")
data("nutrient_lookup_files")


# For example, if we are interested in reaction HMR_4363 in glycolysis pathway:
#HMR_4363: 2-phospho-D-glycerate[c] <=> H2O[c] + PEP[c]
flux_Young_Spleen_HMR_4363<-data.frame(hmr4363=flux_Young_Spleen[grep("HMR_4363",rownames(flux_Young_Spleen)),])
flux_Young_Spleen_HMR_4363







