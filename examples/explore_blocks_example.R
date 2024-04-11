# read data
data("reindeer_ssf")

# spatial stratification
spst <- spat_strat(reindeer, coords = c("x", "y"), colH0 = "original_animal_id")

# explore blocks
explore_blocks(spst)
