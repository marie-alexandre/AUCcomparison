# load required packages ----
if (!require("pacman")) install.packages("pacman") 
pacman::p_load(magrittr, dplyr, usethis, data.table, here)

LOD <- log10(50)

# Creation of HIV_Simu_Dataset_Delta0 ----
HIV_Simu_Dataset_Delta0 <- read.csv2(here::here("data-raw","Individual_dynamics_Delta0.csv"),header = TRUE,stringsAsFactors = FALSE)
usethis::use_data(HIV_Simu_Dataset_Delta0, overwrite = TRUE)

# Creation of HIV_Simu_Dataset_Delta0_cens ----
HIV_Simu_Dataset_Delta0_cens <- HIV_Simu_Dataset_Delta0
HIV_Simu_Dataset_Delta0_cens$cens <- 1*(round(HIV_Simu_Dataset_Delta0_cens$VL,digits=5) <= round(LOD,digits=5))
HIV_Simu_Dataset_Delta0_cens$VL[which(HIV_Simu_Dataset_Delta0_cens$cens == 1)] <- round(LOD,digits=5)
usethis::use_data(HIV_Simu_Dataset_Delta0_cens, overwrite = TRUE)


# Creation of HIV_Simu_Dataset_Delta01 ----
HIV_Simu_Dataset_Delta01 <- read.csv2(here::here("data-raw","Individual_dynamics_Delta01.csv"),header = TRUE,stringsAsFactors = FALSE)
usethis::use_data(HIV_Simu_Dataset_Delta01, overwrite = TRUE)

# Creation of HIV_Simu_Dataset_Delta01_cens ----
HIV_Simu_Dataset_Delta01_cens <- HIV_Simu_Dataset_Delta01
HIV_Simu_Dataset_Delta01_cens$cens <- 1*(round(HIV_Simu_Dataset_Delta01_cens$VL,digits=5) <= round(LOD,digits=5))
HIV_Simu_Dataset_Delta01_cens$VL[which(HIV_Simu_Dataset_Delta01_cens$cens == 1)] <- round(LOD,digits=5)
usethis::use_data(HIV_Simu_Dataset_Delta01_cens, overwrite = TRUE)



# Creation of HIV_Simu_Dataset_Delta025 ----
HIV_Simu_Dataset_Delta025 <- read.csv2(here::here("data-raw","Individual_dynamics_Delta025.csv"),header = TRUE,stringsAsFactors = FALSE)
usethis::use_data(HIV_Simu_Dataset_Delta025, overwrite = TRUE)

# Creation of HIV_Simu_Dataset_Delta025 ----
HIV_Simu_Dataset_Delta025_cens <- HIV_Simu_Dataset_Delta025
HIV_Simu_Dataset_Delta025_cens$cens <- 1*(round(HIV_Simu_Dataset_Delta025_cens$VL,digits=5) <= round(LOD,digits=5))
HIV_Simu_Dataset_Delta025_cens$VL[which(HIV_Simu_Dataset_Delta025_cens$cens == 1)] <- round(LOD,digits=5)
usethis::use_data(HIV_Simu_Dataset_Delta025_cens, overwrite = TRUE)

rm(LOD)