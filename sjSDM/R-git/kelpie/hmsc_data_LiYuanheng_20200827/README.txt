# MODIFY THIS README FILE TO CORRESPOND TO YOUR DATA

# KEEP IT BRIEF; THIS IS NOT INTENDED TO BE A COMPLETE DESCRIPTION OF YOUR DATA

# BUT JUST A BRIEF OVERVIEW TO GET A QUICK IDEA

# FOLLOW THE STRUCTURE AND STYLE OF THIS FILE AS MUCH AS POSSIBLE



This README.text file includes the basic information about the dataset coded as

Li_Yuanheng_Oregon #change this line!



# MAIN QUESTION/AIM OF THE STUDY IN ONE SENTENCE:

To predict the distribution of Malaise-trapped invertebrates in a working forest as a function of remotely-sensed covariates (LiDAR, multispectral).



# SPECIES DATA (Y)

There are two species datasets, and we have made two read_data.R files.

incidence_lidar_mulspec_5_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv
quasiP_lidar_mulspec_5_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200723.csv

The filename beginning with 'incidence' is presence/absence species data
The filename beginning with quasiP is quasi-probability species data [0,1].  We used a DNA spike-in standard to estimate within-species abundances (how a species' biomass changes from sample to sample), and we transformed the read-count data using scales::rescale(Y).

In both datasets there are 1147 arthropod species. The species were extracted from shotgun-sequenced Malaise trapped samples and then subjected to in-silico PCR using the software Kelpie. We used the output to build a reference library of species barcodes (418 bp, BF3BR2 primers), and we mapped the reads from each sample to the reference barcodes. This method largely follows the SPIKEPIPE protocol, except that we used Kelpie to generate the barcode reference dataset (Otso and Nerea are familiar with this paper.)

Because there is no physical PCR and because we set a minimum mapping percent coverage of 50%, false-positive detections should be rare. False negatives, however, are probably pretty common (missing barcodes from the reference and missing detections of known barcodes)



# STUDY DESIGN (S)

There are 94 sampling points (SiteName), and there were two trapping sessions (Session: S1, S2). In some sessions, we used 1 Malaise trap (trap: M1), and in sessions, we used 2 Malaise traps (trap: M1, M2).

The spatial coordinates of each sampling point are Route_x and Route_y

The two trapping sessions caught quite different species because between the two sessions, there was a fire in the general region that caused a lot of smoke in the area and also because of general seasonal turnover.


# COVARIATES (X)

Our primary question is model selection. We have three sets of covariates (GIS, multispectral, LiDAR). Preliminary analyses with NMDS show that species compositions are explained by GIS variables (elevation, canopy.ht, yrs.disturb.min). Our goal is to predict species composition as a function of *only* elevation + multispectral and/or Lidar covariates.

# GIS-type covariates, correlated with landscape structure
elevation (continuous): elevation (meter) of the sampling point
canopy.ht (continuous): mean canopy height (meter) of the sampling point
min.T (continuous): annual minimum temperature
max.T (continuous): annual maximum temperature
precipitation (continuous): annual precipitation
metre.road (continuous): distance to road (meter)
metre.stream (continuous): distance to stream (meter)
yrs.disturb.min (numeric): years since last disturbance
hja (factor): inside or outside the official HJA Andrews Experimental forest border (all sites are in forest)
The above covariates also have scaled version (e.g. elevation.scale)

# multispectral covariates, correlated with vegetation condition and plant species composition
mean.NDVI.scale, mean.EVI.scale, mean.green.scale, mean.bright.scale, mean.wet.scale (continuous):  multispectral satellite data (mean normalized difference vegetation index, enhanced vegetation index, greenness, brightness, wetness), scaled.  The 'mean' refers to the mean value over a area surrounding each sampling point.

# LiDAR covariates, correlated with forest structural complexity
l_Cover_2m_4m.scale, l_Cover_2m_4m_all.scale, l_Cover_2m_max.scale, l_Cover_2m_max_all.scale, l_Cover_4m_16m.scale, l_p25.scale, l_p25_all.scale, l_p95.scale, l_p95_all.scale, l_rumple.scale, l_Cover_2m_4m, l_Cover_2m_4m_all, l_Cover_2m_max, l_Cover_2m_max_all (continuous):  Lidar data, scaled.


# unscaled versions of LiDAR and multispectral covariates
l_Cover_4m_16m, l_p25, l_p25_all, l_p95, l_p95_all, l_rumple:  Lidar data, unscaled
mean.NDVI, mean.EVI, mean.bright, mean.green, mean.wet:  multispectral satellite data, unscaled


# HOW TO SET UP A REASONABLE PILOT MODEL FOR XFormula:

Limit dataset to M1,S1 (Malaise trap 1 in Session 1) for training. M2,S1 (Malaise trap 2 in Session 1) can be validation dataset for Session 1.

Pilot model 1 (GIS):  Route_x, Route_y, elevation, canopy.ht, yrs.disturb.min
Pilot model 2 (multispectral):  Route_x, Route_y, mean.NDVI.scale, mean.EVI.scale, mean.green.scale, mean.bright.scale, mean.wet.scale
Pilot model 3 (LiDAR):  Route_x, Route_y, l_Cover_2m_4m.scale, l_Cover_4m_16m.scale, l_p25.scale, l_p95.scale, l_rumple.scale


Transformations suggested for continuous variables:  no suggestions


# TRAITS (Tr)

NO



# HOW TO SET UP A REASONABLE PILOT MODEL FOR TrFormula:

NO



# PHYLOGENY (P)

NO



# CAN WE USE YOUR DATA AS EXAMPLE IN THE COURSE: YES # change to NO as needed (?)



# DO YOU POTENTIALLY WISH TO COLLABORATE WITH YOUR DATA ANALYSES AT THE LEVEL OF CO-AUTHORSHIP: YES # change to "NO" or "MAYBE" as needed (?)

The default plan is that we will advise you during the course with your data analyses, including setting up and fitting the models, interpreting the results, and writing the paper. If you find our advice useful, we would like you to mention that in the acknowledgements.
If you wish, we may also be able to take a more active role in your manuscript preparation, so that a co-authorship becomes justified. If you answer "YES", this means that in principle you would welcome such a contribution. Of course, whether an instructor from the course eventually becomes a co-author or not will depend on their actual level of contribution to the manuscript.
