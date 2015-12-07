library(ape)
library(rncl)
library(rotl)
library(phylobase)
library(phytools)
library(rgbif)
library(ggplot2)
library(maps)

### Basic plotting -------------------------------------------------------------

## labelling nodes according to their support values

tr <- rncl::read_nexus_phylo("data/20091017.Antarctic.nex.con")[[1]]

plot(ladderize(tr), show.tip.label=FALSE)
nodelabels(text=tr$node.label, frame="rect", cex=.5, bg="white")

node_colors <- sapply(tr$node.label, function(x) {
    if(is.na(x)) return(NA)
    if(x < .75) return("white")
    if(x >= .75 && x < .95) return("yellow")
    if(x >= .95) return("red")
    })


plot(ladderize(tr), show.tip.label=FALSE)
nodelabels(text=rep("", tr$Nnode), frame="circ", bg=node_colors, cex=.3)

## Zooming on one part of the tree

plot(ladderize(tr), show.tip.label=FALSE)
zoom(ladderize(tr), grep("moseleyi", tr$tip.label), subtree=TRUE)


### Tree and data --------------------------------------------------------------

## `phylobase` for associating data to tips + phylogenetic

data(geospiza_raw)
geo_ape <- geospiza_raw$tree
geo_data <- geospiza_raw$data

geo <- phylo4d(geo_ape, geo_data,missing="OK")
geo <- prune(geo, "olivacea")

## independent constrasts

geo_tr <- as(extractTree(geo), "phylo")

beak <- tdata(geo, "tip")$beakD
names(beak) <- tipLabels(geo)
wing <- tdata(geo, "tip")$wingL
names(wing) <- tipLabels(geo)

par(mfrow = c(1, 2))
contMap(geo_tr, wing)
contMap(geo_tr, beak)

plot(beakD ~ wingL, data = geo_data)
summary(lm(beakD ~ wingL - 1, data = geo_data))

pic_beak <- pic(beak, geo_tr)
pic_wing <- pic(wing, geo_tr)

plot(pic_wing, pic_beak)
summary(lm(pic_beak ~ pic_wing - 1))

par(mfrow = c(1, 2))
plot(beakD ~ wingL, data = geo_data, main = "raw")
plot(pic_wing, pic_beak, pch = 15, main = "pic")

## Simulation of continuous traits under Brownian motion

tr <- as(rcoal(12), "phylo4")
vmat <- as(tr, "phylo4vcov")
mat <- cov2cor(vmat)
library(MASS)
trvec <- mvrnorm(5, mu=rep(0,12), Sigma=vmat)
treed <- phylo4d(tr, tip.data=t(data.frame(trvec)))
plot(treed)

### Getting data from GenBank --------------------------------------------------

gb <- read.csv(file="data/OLoughlinAppendixTable.csv", stringsAsFactors=FALSE)
mb <- read.GenBank(gb$MB[nzchar(gb$MB)])
write.dna(mb, file="mb.fas", format="fasta", nbcol=1, colw=50)
system("muscle -in mb.fas -out mb.afa")


### rotl (Open Tree of Life ) --------------------------------------------------

## you might need the development version of the package for this
## demonstration to work install devtools and then:
## devtools::install_github("ropensci/rotl")

taxa <- tnrs_match_names(names = c("Escherichia colli",
                                   "Chlamydomonas reinhardtii",
                                   "Drosophila melanogaster",
                                   "Arabidopsis thaliana",
                                   "Rattus norvegicus",
                                   "Mus musculus",
                                   "Cavia porcellus",
                                   "Xenopus laevis",
                                   "Saccharomyces cervisae",
                                   "Danio rerio"))

tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]])
plot(tree, cex = .8, label.offset = .1, no.margin = TRUE)

cat_studies <- studies_find_studies(property = "ot:focalCladeOTTTaxonName",
                                    value = "felidae")
cat_studies

cat_tree <- get_study_tree(study_id = cat_studies[["study_ids"]][1],
                           tree_id = cat_studies[["tree_ids"]][1])
cat_tree

plot(cat_tree, label.offset=.0001, cex=.75, no.margin=TRUE)


## Import the felid tree using study and tree IDs discovered with
## studies_find_studies() in the manuscript
cat_tree <- get_study_tree(study_id ="pg_1981",
                           tree_id = "tree4052")

## Find the species of Lynx in the phylogeny
cat_species <- cat_tree$tip.label
lynx_species <- grep("^Lynx", cat_tree$tip.label, value = TRUE)

## Match the Lynx species to the GBIF identifiers
gbif_keys <- sapply(lynx_species,
                    function(x) name_backbone(name = x)$speciesKey,
                    USE.NAMES = FALSE)

## Search for the GBIF records for these species
lynx_loc <- occ_search(taxonKey = gbif_keys, limit = 500,
                       return = "data", fields = "minimal",
                       hasCoordinate = TRUE)

## Make a data frame of the results
lynx_loc <- do.call("rbind", lynx_loc)
names(lynx_loc)[1] <- "Species"

## Clean up the data with missing locality data
lynx_loc[["decimalLatitude"]] <- as.numeric(lynx_loc[["decimalLatitude"]])
lynx_loc[["decimalLongitude"]] <- as.numeric(lynx_loc[["decimalLongitude"]])
lynx_loc[lynx_loc[["decimalLatitude"]] == 0 &
         lynx_loc[["decimalLongitude"]] == 0,
         c("decimalLatitude","decimalLongitude")] <- c(NA, NA)
lynx_loc <- lynx_loc[complete.cases(lynx_loc), ]

## Draw the map
world <- map_data("world")

ggplot(lynx_loc) +
    annotation_map(world, fill="gray40", color="gray40") +
    geom_point(aes(y = decimalLatitude, x = decimalLongitude, color = Species),
               size = 1) +
    coord_map(projection = "mercator", orientation = c(90, 0, 0)) +
    xlab("Longitude") + ylab("Latitude") +
    theme(legend.position="top", legend.key = element_rect(fill = "gray40")) +
    ylim(c(0,72))
